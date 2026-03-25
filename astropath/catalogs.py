"""
Methods related to public catalogs for PATH analysis

This module provides functionality to query public astronomical surveys
and prepare catalogs for PATH analysis.
"""
import os

import numpy as np

from astropy.coordinates import SkyCoord
from astropy import units
from astropy.table import Table

from astroquery.vizier import Vizier

from IPython import embed


def query_catalog(survey: str, coord: SkyCoord, ssize: float,
                  debug: bool = False, skip_NGC: bool = False,
                  dust_correct: bool = True,
                  star_galaxy_sep: dict = None):
    """Query a public catalog for PATH analysis.

    The catalog is cleaned up to remove objects that are too bright
    or likely to be stars.

    Args:
        survey (str): Name of the survey ['Pan-STARRS', 'DECaL']
        coord (astropy.coordinates.SkyCoord): Coordinates of the transient
        ssize (float): Radius of the query in arcmin
        debug (bool, optional): Show some plots. Defaults to False.
        skip_NGC (bool, optional): Skip adding NGC galaxies. Defaults to False.
        dust_correct (bool, optional): Dust correct the magnitudes. Defaults to True.
        star_galaxy_sep (dict, optional): Additional star/galaxy separation parameters.
            For Pan-STARRS, can include 'PS1_PSC_cut' (default 0.83)
            and 'PS1_PSC_func' (callable to get PSC values).

    Raises:
        ValueError: If survey is not supported

    Returns:
        tuple: 3 items --
            (astropy.table.Table) table of galaxies,
            (str) Magnitude key (e.g. 'Pan-STARRS_r'),
            (astropy.table.Table) table of stars rejected
    """
    # Import here to allow astropath to work without frb package
    try:
        from frb.surveys import survey_utils
    except ImportError:
        raise ImportError(
            "The 'frb' package is required for catalog queries. "
            "Install with: pip install frb"
        )

    # Survey specific queries
    if survey == 'Pan-STARRS':
        query_fields = ['rPSFLikelihood', 'iPSFMag', 'iKronMag']
    elif survey == 'DECaL':
        query_fields = ['shape_r', 'type']
    else:
        query_fields = None
    print(f"Using: {survey} survey")

    survey_obj = survey_utils.load_survey_by_name(
        survey, coord, ssize * units.arcmin)

    # Grab the catalog
    catalog = survey_obj.get_catalog(query_fields=query_fields)

    # Empty?
    if len(catalog) == 0:
        print(f"No objects in the catalog of your survey within {ssize} arcmin")
        return catalog, None, None

    # Clean up the catalog
    if survey == 'Pan-STARRS':
        catalog, mag_key, stars = _clean_panstarrs_catalog(
            catalog, coord, ssize, star_galaxy_sep)
    elif survey == 'DECaL':
        catalog, mag_key, stars = _clean_decal_catalog(catalog)
    else:
        raise ValueError(f"Not ready for this survey: {survey}")

    # Dust correct?
    if dust_correct and len(catalog) > 0:
        catalog = _apply_dust_correction(catalog, mag_key)

    # Add nearby NGC galaxies
    if not skip_NGC and len(catalog) > 0:
        add_NGC_galaxies(coord, catalog, mag_key)

    if debug:
        _debug_plots(catalog, survey)

    # Return
    return catalog, mag_key, stars


def _clean_panstarrs_catalog(catalog: Table, coord: SkyCoord,
                              ssize: float, star_galaxy_sep: dict = None):
    """Clean Pan-STARRS catalog for PATH analysis.

    Args:
        catalog: Raw Pan-STARRS catalog
        coord: Central coordinates
        ssize: Search radius in arcmin
        star_galaxy_sep: Dict with optional 'PS1_PSC_cut' and 'PS1_PSC_func'

    Returns:
        tuple: (cleaned catalog, mag_key, stars table)
    """
    # Default star/galaxy separation parameters
    if star_galaxy_sep is None:
        star_galaxy_sep = {}
    PS1_PSC_cut = star_galaxy_sep.get('PS1_PSC_cut', 0.83)
    PS1_PSC_func = star_galaxy_sep.get('PS1_PSC_func', None)

    # Remove non-detections in the r-band
    not_junk = catalog['Pan-STARRS_r'] > 0.
    catalog = catalog[not_junk]

    if len(catalog) == 0:
        return catalog, 'Pan-STARRS_r', Table()

    # Various star/galaxy cuts
    cut_size = catalog['rKronRad'] > 0.
    cut_mag = catalog['Pan-STARRS_r'] > 15.
    cut_point = np.log10(np.abs(catalog['rPSFLikelihood'])) < (-2)
    cut_psf = (catalog['iPSFMag'] - catalog['iKronMag_1']) > 0.05

    # PS1 PSC cut (if function provided)
    if PS1_PSC_func is not None:
        PS1_PSC_func(catalog)
        cut_psc = catalog['PS1_PSC'] < PS1_PSC_cut
    else:
        # Without PS1_PSC data, skip this cut
        cut_psc = np.ones(len(catalog), dtype=bool)

    # Cross-match to remove GAIA stars with Vizier
    cut_gaia_stars, gaia_catalog = _remove_gaia_stars(catalog, coord, ssize)

    keep = cut_size & cut_mag & cut_point & cut_psc & cut_psf & cut_gaia_stars

    # Convert all "stars" fainter than mag = 20 to a galaxy
    too_faint = (catalog['Pan-STARRS_r'] > 20.) & np.isfinite(catalog['Pan-STARRS_r'])
    keep[too_faint] = True

    # Half-light radius
    mag_key = 'Pan-STARRS_r'
    catalog['ang_size'] = catalog['rKronRad'].copy()
    catalog['ID'] = catalog['Pan-STARRS_ID'].copy()

    # Cut
    stars = catalog[~keep]
    catalog = catalog[keep]

    return catalog, mag_key, stars


def _clean_decal_catalog(catalog: Table):
    """Clean DECaL catalog for PATH analysis.

    Args:
        catalog: Raw DECaL catalog

    Returns:
        tuple: (cleaned catalog, mag_key, stars table)
    """
    mag_key = 'DECaL_r'

    # Cuts
    cut_mag = (catalog[mag_key] > 14.) & np.isfinite(catalog[mag_key])
    cut_more_stars = catalog['DECaL_type'] != 'PSF'

    keep = cut_mag & cut_more_stars

    # Take everything fainter than 23 mag as a galaxy
    too_faint = (catalog[mag_key] > 23.) & np.isfinite(catalog[mag_key])
    keep[too_faint] = True

    # Half-light radius
    catalog['ang_size'] = catalog['shape_r']
    bad_ang = catalog['ang_size'] == 0.
    if np.any(bad_ang) > 0:
        print(f"WARNING: Found {np.sum(bad_ang)} objects with zero ang_size. Setting to 0.7 arcsec")
    catalog['ang_size'][bad_ang] = 0.7

    catalog['ID'] = catalog['DECaL_ID'].copy()

    # Cut
    stars = catalog[~keep]
    catalog = catalog[keep]

    return catalog, mag_key, stars


def _remove_gaia_stars(catalog: Table, coord: SkyCoord, ssize: float):
    """Cross-match with Gaia DR3 to identify stars.

    Args:
        catalog: Catalog with 'Pan-STARRS_ID' column
        coord: Central coordinates
        ssize: Search radius in arcmin

    Returns:
        np.ndarray: Boolean mask where True = not a Gaia star
        astropy.table.Table: Gaia catalog
    """
    Vizier.ROW_LIMIT = -1
    query_radius = ssize / 60. * units.deg

    try:
        vizier = Vizier(columns=["*", "+PSS", "+PGal"], catalog='I/355/gaiadr3')
        gaia_catalog = vizier.query_region(coord, radius=query_radius)[0]#, catalog='I/355/gaiadr3')

        # Cut on PSS
        gaia_catalog = gaia_catalog[gaia_catalog['PSS'] > 0.99]

        #gaia_catalog = Vizier.query_region(
        #    coord, radius=query_radius, catalog='I/355/gaiadr3')
        if len(gaia_catalog) > 0:
            gaia_star_ids = gaia_catalog[0]['PS1']  # Pan-STARRS IDs matched to Gaia stars
            cut_gaia_stars = np.logical_not(
                np.isin(catalog['Pan-STARRS_ID'], gaia_star_ids))
        else:
            cut_gaia_stars = np.ones(len(catalog), dtype=bool)
    except Exception as e:
        print("Vizier Gaia query failed:", e)
        cut_gaia_stars = np.ones(len(catalog), dtype=bool)

    return cut_gaia_stars, gaia_catalog


def _apply_dust_correction(catalog: Table, mag_key: str):
    """Apply dust correction to magnitudes.

    Args:
        catalog: Catalog with magnitude column
        mag_key: Name of magnitude column

    Returns:
        Table: Catalog with dust-corrected magnitudes
    """
    try:
        from frb.galaxies import nebular
        from frb.galaxies import photom as frbphotom
    except ImportError:
        print("WARNING: frb package not available for dust correction. Skipping.")
        return catalog

    # Get E(B-V) at first source position
    coord = SkyCoord(ra=catalog['ra'][0],
                     dec=catalog['dec'][0], unit='deg',
                     frame='icrs')
    EBV = nebular.get_ebv(coord)['meanValue']
    dust_correct = frbphotom.extinction_correction(
        mag_key, EBV, required=True)
    mag_dust = 2.5 * np.log10(1. / dust_correct)
    catalog[mag_key] += mag_dust

    return catalog


def add_NGC_galaxies(coord: SkyCoord, catalog: Table, mag_key: str,
                     sep_arcmins: float = 60):
    """Add NGC galaxies to the table of candidates.

    The input table is modified in place.

    Args:
        coord (SkyCoord): coordinates of the FRB
        catalog (Table): table of existing candidates
        mag_key (str): name of the magnitude key
        sep_arcmins (float, optional): max radius to search to, arcmins
    """
    try:
        from pyongc.ongc import nearby
    except ImportError:
        print("WARNING: pyongc not available. Skipping NGC galaxy addition.")
        return

    print("Add bright NGC/IC galaxies")
    ra_str = coord.ra.to_string(units.hour, sep=':', pad=True, precision=2)
    dec_str = coord.dec.to_string(units.degree, alwayssign=True, sep=':', pad=True, precision=2)
    print("FRB Position: {}, {}".format(ra_str, dec_str))

    nearby_ngcs = nearby('{} {}'.format(ra_str, dec_str), separation=sep_arcmins)
    if len(nearby_ngcs) == 0:
        print("No NGC/IC galaxies within {} arcmins of the FRB".format(sep_arcmins))

    for ngc_tuple in nearby_ngcs:
        separation = ngc_tuple[1] * 60.  # arcmin

        ngc = ngc_tuple[0]
        print("{0} is {1:.2f} arcmin from the FRB".format(ngc.name, separation))
        coord_ngc = SkyCoord(ngc.ra, ngc.dec, frame='icrs',
                             unit=(units.hour, units.deg))
        ra = coord_ngc.ra.deg
        dec = coord_ngc.dec.deg

        # ngc.dimensions returns: (MajAx, MinAx, P.A.) in arcmins
        dimensions = ngc.dimensions
        if dimensions[0] is None or dimensions[1] is None:
            continue
        else:
            ang_size = np.max([dimensions[0], dimensions[1]]) * 60.  # arcsec

        # ngc.magnitudes returns: (Bmag, Vmag, Jmag, Hmag, Kmag)
        # Vmag is 500-700 nm, close to Pan-STARRS r-mag centered at ~621 nm
        magnitudes = np.array(ngc.magnitudes)
        vmag = magnitudes[1]
        if vmag is None:
            print("Vmag is None for {}, trying to find another magnitude".format(ngc.name))
            not_none = magnitudes != None
            if np.any(not_none):
                vmag = np.nanmean(magnitudes[not_none])
            else:
                print("No magnitudes available for {}, skipping".format(ngc.name))
                continue

        # Add the above values to the catalog
        catalog.add_row()
        catalog['ra'][-1] = ra
        catalog['dec'][-1] = dec
        catalog['ang_size'][-1] = ang_size
        catalog[mag_key][-1] = vmag


def _debug_plots(catalog: Table, survey: str):
    """Show debug plots for catalog."""
    try:
        import seaborn as sns
        from matplotlib import pyplot as plt
    except ImportError:
        print("seaborn/matplotlib not available for debug plots")
        return

    if survey == 'Pan-STARRS':
        sns.histplot(x=catalog['Pan-STARRS_r'])
        plt.show()
        sns.histplot(x=catalog['rKronRad'])
        plt.show()
        sns.histplot(x=catalog['rPSFLikelihood'])
        plt.show()
