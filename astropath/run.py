"""
Run PATH on a single FRB given a dictionary of inputs

This module provides a high-level interface for running PATH analysis.
"""

import numpy as np

from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy import units

from astropath import catalogs
from astropath import path

from IPython import embed
import time


def run_on_dict(idict: dict,
                verbose: bool = False,
                verbose_sims: bool = False,
                catalog: Table = None,
                mag_key: str = None,
                skip_NGC: bool = False,
                dust_correct: bool = True,
                use_local: bool = False,
                star_galaxy_sep: dict = None,
                sep_cull: bool = False):
    """Run PATH on a single FRB, given a dictionary of inputs.

    The idict must have the following keys:
        - ra (float): Right ascension (deg)
        - dec (float): Declination (deg)
        - ssize (float): Radius for probing the survey (arcmin)
            This should usually be set with the set_anly_sizes() method
        - ltype (str): Localization type (e.g. 'eellipse', 'healpix')
        - lparam (dict): Localization parameters
            For 'eellipse': {'a': 0.2, 'b': 0.2, 'theta': 0.}
                * a, b are in units of arcsec
                * theta is in deg and is E from N
            For 'healpix': {'healpix_data': healpix_data, 'healpix_nside': nside,
                            'healpix_coord': coord_sys, 'healpix_ordering': ordering}
        - priors (dict): PATH priors
            e.g. {'P_O_method': 'inverse', 'PU': 0.5, 'scale': 0.5,
                  'theta_PDF': 'exp', 'theta_max': 6.}
            Can also include 'survey' key for automatic catalog query.
        - max_box (float): Maximum box size for PATH analysis (arcsec)
            This should usually be set with the set_anly_sizes() method

    The dict needs to be simple enough to serialize with JSON for eellipse localization.

    Args:
        idict (dict): Input dictionary with FRB and analysis parameters
        verbose (bool, optional): Print more diagnostic info. Defaults to False.
        catalog (astropy.table.Table, optional): Catalog of candidates.
            If None and idict['priors']['survey'] is set, will query the survey.
            Must have columns: 'ra', 'dec', 'ang_size', and the magnitude column
            specified by mag_key. Optionally includes 'ID' for source identification.
        mag_key (str, optional): Magnitude column key for the catalog.
            If None and catalog is queried, will be determined from the survey.
        skip_NGC (bool, optional): Skip adding NGC galaxies when querying. Defaults to False.
        dust_correct (bool, optional): Dust correct magnitudes when querying. Defaults to True.
        star_galaxy_sep (dict, optional): Star/galaxy separation parameters for catalog query.
        use_local (bool, optional): Use the local method for calculating posteriors. Defaults to False.
            If True, will use the local method for calculating posteriors.
            If False, will use the fixed method for calculating posteriors.

    Returns:
        tuple: 6 items --
            - pandas.DataFrame: Candidates with computed posteriors
            - float: P(U|x) - probability of unseen host
            - PATH: PATH object used for analysis
            - str: mag_key used
            - astropy.Table: Catalog cut down to the analysis box
            - astropy.Table or None: Stars rejected from catalog (if queried), else None
    """
    # Unpack coordinates
    coord = SkyCoord(ra=idict['ra'], dec=idict['dec'], unit='deg')

    # Query catalog if not provided and survey is specified
    stars = None
    if catalog is None:
        survey = idict.get('priors', {}).get('survey')
        if survey is not None:
            catalog, mag_key, stars = catalogs.query_catalog(
                survey,
                coord,
                idict['ssize'],
                skip_NGC=skip_NGC,
                dust_correct=dust_correct,
                star_galaxy_sep=star_galaxy_sep
            )
        else:
            raise ValueError(
                "catalog is required. Pass an astropy.table.Table with "
                "columns: 'ra', 'dec', 'ang_size', and the magnitude column, "
                "OR set idict['priors']['survey'] to query automatically."
            )

    # Validate mag_key
    if mag_key is None:
        raise ValueError("mag_key is required. Specify the column name for magnitudes.")

    # Localization
    if idict['ltype'] == 'eellipse':
        eellipse = idict['lparam']
    elif idict['ltype'] == 'healpix':
        # For healpix, lparam should contain healpix_data and related parameters
        pass
    else:
        raise ValueError(f"Unsupported localization type: {idict['ltype']}. "
                        f"Supported: 'eellipse', 'healpix'")

    # Handle empty catalog
    if len(catalog) == 0:
        return None, None, None, None, None, None
    
    # Set boxsize according to the largest galaxy (arcsec)
    #  or the survey size, whichever is larger
    # TODO -- X decided this is not worth the extra compute
    #cat_coord = SkyCoord(ra=catalog['ra'].data, dec=catalog['dec'].data, unit='deg')
    #sep = coord.separation(cat_coord).to('arcsec')
    #max_sep = np.max(sep + 6*catalog['ang_size'].data*units.arcsec)
    #embed(header='run.py:117')

    #box_hwidth = max(idict['max_box'],  # This should be the survey size
    #    max_sep.value)
    box_hwidth = idict['max_box']  # This should be the survey size

    # Calculate separation from localization center
    catalog_coord = SkyCoord(ra=catalog['ra'].data, dec=catalog['dec'].data, unit='deg')
    sep = coord.separation(catalog_coord).to('arcsec')
    catalog['sep'] = sep.arcsec
    # Cut down the catalog based on ssize (usually this should do nothing)
    keep = sep < idict['ssize']*60*units.arcsec
    cut_catalog = catalog[keep]

    if len(cut_catalog) == 0:
        return None, None, None, None, None, None

    # Initialize PATH
    Path = path.PATH()

    # Set up localization
    if idict['ltype'] == 'eellipse':
        Path.init_localization('eellipse',
                               center_coord=coord,
                               eellipse=eellipse)
    elif idict['ltype'] == 'healpix':
        lparam = idict['lparam']
        Path.init_localization('healpix',
                               center_coord=coord,
                               healpix_data=lparam.get('healpix_data'),
                               healpix_nside=lparam.get('healpix_nside'),
                               healpix_coord=lparam.get('healpix_coord', 'C'),
                               healpix_ordering=lparam.get('healpix_ordering', 'RING'))

    # Initialize candidates
    Path.init_candidates(cut_catalog['ra'],
                         cut_catalog['dec'],
                         cut_catalog['ang_size'],
                         mag=cut_catalog[mag_key],
                        )

    # Add ID if available
    if 'ID' in cut_catalog.colnames:
        Path.candidates['ID'] = cut_catalog['ID']

    # Add separation from localization center
    ccand = SkyCoord(ra=Path.candidates['ra'], dec=Path.candidates['dec'], unit='deg')
    sep = ccand.separation(coord)
    Path.candidates['sep'] = sep.arcsec

    # Set up priors from idict
    priors_dict = idict['priors']

    # Candidate prior
    P_O_method = priors_dict.get('P_O_method', 'inverse')
    P_U = priors_dict['PU']
    Path.init_cand_prior(P_O_method, P_U=P_U)

    # Offset prior
    theta_PDF = priors_dict['theta_PDF']
    theta_max = priors_dict.get('theta_max', 6.)
    scale = priors_dict['scale']
    Path.init_theta_prior(theta_PDF, theta_max, scale)

    # Calculate priors
    p_O = Path.calc_priors()

    # Calculate step size based on localization and galaxy sizes
    if idict['ltype'] == 'eellipse':
        if idict['pmode'] == 'fixed':
            min_ang = np.nanmin(cut_catalog['ang_size'].data)
            if eellipse['b'] > min_ang:
                correction = 'p_wO'
                step_size = eellipse['b'] / 20.
            else:
                correction = 'L_wx'
                step_size = min_ang / 20.
            idict['step_size'] = step_size
        elif idict['pmode'] == 'local':
            assert 'step_size' in idict, "step_size is required for local mode"
            assert 'step_size_mode' in idict, "step_size_mode is required for local mode"
            correction = None
    else:
        # For healpix, use default step size
        raise ValueError("Healpix localization is not supported yet.")

    # Add to dict
    if 'index' in idict:
        index = idict['index']

    # Calculate posteriors
    start_time = time.perf_counter()
    if idict['pmode'] == 'local':
        if verbose_sims and ('index' in idict):
            max_ang = np.nanmax(cut_catalog['ang_size'].data)
            step_size = idict['step_size']
            print(f'FRB {index}: local, ssize={box_hwidth}, max_box_hwidth={max_ang*6.}, step_size={step_size}')
        P_Ox, P_Ux = Path.calc_posteriors(
            'local',
            survey_radius=idict['ssize']*60,
            step_size=idict['step_size'],
            step_size_mode=idict['step_size_mode'],
            sep_cull=sep_cull)
    elif idict['pmode'] == 'fixed':
        # Memory check
        if idict['ssize']*60 / step_size > 10000:
            raise ValueError(f"Fixed grid would be {int(idict['ssize']*60 / step_size)} pixels.  \nYour array will need >100Gb RAM.  Try local or reduce your ssize if you can")
        if verbose_sims and ('index' in idict):
            step_size = idict['step_size']
            max_ang = np.nanmax(cut_catalog['ang_size'].data)
            print(f'FRB {index}: fixed, correction: {correction}, ssize={box_hwidth}, max_box_hwidth={max_ang*6.}, step_size={step_size}')
        P_Ox, P_Ux = Path.calc_posteriors('fixed',
                                       box_hwidth=box_hwidth,
                                       survey_radius=idict['ssize']*60,
                                       step_size=idict['step_size'],
                                       use_numba=idict['use_numba'],
                                       correction=correction, 
                                       verbose_sims=verbose_sims)
    else:
        raise ValueError(f"Unsupported posterior mode: {idict['pmode']}. "
                        f"Supported: 'local', 'fixed'")

    end_time = time.perf_counter()
    if verbose_sims and ('index' in idict):
        print(f'FRB {index}: execution time: {end_time - start_time:.4f} seconds')

    #embed(header='run.py:231')

    # Add photo-z columns if available in catalog
    photoz_columns = ['z_phot_median', 'z_phot_l68', 'z_phot_u68',
                      'z_phot_l95', 'z_phot_u95', 'z_spec',
                      'z_phot', 'z_photErr']
    for key in photoz_columns:
        if key in cut_catalog.colnames:
            Path.candidates[key] = cut_catalog[key]

    # Sort by posterior probability
    Path.candidates.sort_values(by='P_Ox', ascending=False, inplace=True)
    top_cand_POx = Path.candidates.iloc[0]['P_Ox']
    if verbose_sims and ('index' in idict):
        print(f'FRB {index}: P(O1|x)={top_cand_POx:.5f}, P(U|x)={P_Ux:.5f}')

    # Print results if verbose
    if verbose:
        print_cols = ['ra', 'dec', 'ang_size', 'mag', 'P_O', 'P_Ox']
        if 'ID' in Path.candidates.columns:
            print_cols = ['ID'] + print_cols
        print(Path.candidates[print_cols])
        print(f"P_Ux = {P_Ux}")

    # Return
    return Path.candidates, P_Ux, Path, mag_key, cut_catalog, stars


def set_anly_sizes(ltype: str, lparam: dict):
    """Set the analysis sizes for PATH.

    This helper function computes appropriate survey search radius (ssize)
    and maximum analysis box size (max_box) based on localization parameters.

    Note: These are guidelines recommended by the PATH
    developers.  You should use your own judgement

    Args:
        ltype (str): Type of localization ['eellipse', 'healpix']
        lparam (dict): Parameters for localization
            For 'eellipse': {'a': semi_major_arcsec, 'b': semi_minor_arcsec, 'theta': PA_deg}
            For 'healpix': Should compute from healpix resolution

    Returns:
        tuple: (ssize, max_box)
            ssize (float): Radius of box to probe the survey (arcmin)
            max_box (float): Maximum box size for PATH analysis (arcsec)

    Raises:
        ValueError: If ltype is not supported
    """
    if ltype == 'eellipse':
        # Survey search size based on semi-major axis
        ssize = 10 * lparam['a'] / 60.  # arcmin
    elif ltype == 'healpix':
        # For healpix, estimate from nside if available
        if 'healpix_nside' in lparam:
            # Pixel size in arcmin for given nside
            nside = lparam['healpix_nside']
            pixel_size_arcmin = 60. * np.degrees(np.sqrt(4 * np.pi / (12 * nside**2)))
            ssize = 10 * pixel_size_arcmin
        else:
            # Default fallback
            ssize = 5.0  # arcmin
    else:
        raise ValueError(f"Unsupported localization type: {ltype}. "
                        f"Supported: 'eellipse', 'healpix'")

    # Enforce minimum search radius
    ssize = max(ssize, 10.) # No less than 10 arcmin

    # Max box for PATH analysis (in arcsec)
    max_box = ssize * 60.  # arcsec

    return ssize, max_box


def build_idict(ra: float, dec: float,
                ltype: str = 'eellipse',
                lparam: dict = None,
                P_O_method: str = 'inverse',
                PU: float = 0.0,
                scale: float = 0.5,
                theta_PDF: str = 'exp',
                theta_max: float = 6.0,
                use_numba: bool = False,
                pmode: str = 'local',
                step_size_mode: str = 'relative',
                step_size: float = 0.05,
                survey: str = None,
                ssize: float = None,
                max_box: float = None):
    """Build an input dictionary for run_on_dict().

    This is a convenience function to construct the input dictionary
    with proper structure and optional auto-calculation of analysis sizes.

    Args:
        ra (float): Right ascension in degrees
        dec (float): Declination in degrees
        ltype (str): Localization type. Default 'eellipse'
        lparam (dict): Localization parameters.
            For 'eellipse': {'a': arcsec, 'b': arcsec, 'theta': deg}
            Required if ltype is 'eellipse'.
        P_O_method (str): Method for candidate prior. Default 'inverse'
        PU (float): Prior probability of unseen host. Default 0.0
        scale (float): Scale factor for offset prior. Default 0.5
        theta_PDF (str): PDF for offset prior. Default 'exp'
        theta_max (float): Maximum offset for prior. Default 6.0
        survey (str, optional): Survey name for automatic catalog query.
            Supported: 'Pan-STARRS', 'DECaL'. If None, catalog must be provided
            to run_on_dict().
        use_numba (bool): Use numba for calculations. Default False.
        ssize (float): Survey search radius in arcmin. If None, auto-computed.
        max_box (float): Maximum analysis box in arcsec. If None, auto-computed.

    Returns:
        dict: Input dictionary suitable for run_on_dict()

    Raises:
        ValueError: If required parameters are missing
    """
    if lparam is None:
        raise ValueError("lparam is required. For 'eellipse', provide "
                        "{'a': arcsec, 'b': arcsec, 'theta': deg}")

    # Auto-compute analysis sizes if not provided
    if ssize is None or max_box is None:
        auto_ssize, auto_max_box = set_anly_sizes(ltype, lparam)
        if ssize is None:
            ssize = auto_ssize
        if max_box is None:
            max_box = auto_max_box

    priors = {
        'P_O_method': P_O_method,
        'PU': PU,
        'scale': scale,
        'theta_PDF': theta_PDF,
        'theta_max': theta_max,
    }
    if survey is not None:
        priors['survey'] = survey

    idict = {
        'ra': ra,
        'dec': dec,
        'ssize': ssize, # arcmin
        'ltype': ltype,
        'lparam': lparam,
        'use_numba': use_numba,
        'pmode': pmode,
        'step_size_mode': step_size_mode,
        'step_size': step_size,
        'priors': priors,
        'max_box': max_box, # arcsec
    }

    return idict
