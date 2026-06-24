
import numpy as np
import pandas
import multiprocessing

from astropy.coordinates import SkyCoord
from astropy import units
from astropy.coordinates import search_around_sky
from astropy.table import Table

from astropath.run import run_on_dict, set_anly_sizes

from IPython import embed

import gc

    

def run_dict_wrapper(args):
    """
    Run the simulation on a dictionary of parameters.

    Mainly used for multiprocessing.

    Args:
        idx (int): The index of the simulation.
        idict (dict): The dictionary of simulation parameters.
        catalog (pandas.DataFrame): The catalog of data.

    Returns:
        pandas.DataFrame: The sorted table of candidates.
        int: The index of the simulationl for book-keeping.
        Path: The Path object
    """
    idx, idict, catalog = args
    catalog = Table.from_pandas(catalog)
    candidates, P_Ux, Path, mag_key, cut_catalog, stars = \
        run_on_dict(idict, catalog=catalog, mag_key='mag', sep_cull=True, verbose_sims=False)

    if candidates is None or len(candidates) == 0:
        return idx, None

    sv_tbl = candidates.sort_values('P_Ox', ascending=False)
    sv_tbl['gal_ID'] = sv_tbl.index.values
    # Return ONLY what full() consumes; no Path, no object/SkyCoord columns
    keep = [c for c in ['ra', 'dec', 'ang_size', 'mag', 'ID', 'sep',
                        'P_O', 'P_Ox', 'P_Ux', 'gal_ID'] if c in sv_tbl.columns]
    return idx, sv_tbl[keep].copy()
    

def full(frbs:pandas.DataFrame, catalog:pandas.DataFrame, 
             prior_dict:dict,
             multi:bool=True,
             ncpu:int=4,
             debug:bool=False):
    """
    Run the PATH simulation with a given catalog and priors.

    This function will generate the FRB dicts, build the galaxy catalog tables,
    and run the PATH simulation.
    The simulation results will be returned as a pandas dataframe.

    The method uses multiprocessing to run the PATH simulation on multiple FRBs in parallel
    if multi is True.

    The catalog dataframe must have the following columns:
    - ra (float): The right ascension of the galaxy
    - dec (float): The declination of the galaxy
    - ang_size (float): The angular size of the galaxy
    - mag (float): The magnitude of the galaxy
    - ID (int): The ID of the galaxy

    The prior dictionary must have the following keys:
    - P_O_method (str): The method to use for the prior on the host
    - PU (float): The prior on the unseen host
    - scale (float): The scale of the prior
    - theta_PDF (str): The PDF to use for the prior on the theta
    - theta_max (float): The maximum value of the theta prior

    Args:
        frbs(pandas.DataFrame): The FRBs dataframe
        catalog(pandas.DataFrame): The catalog dataframe
        prior_dict (dict): The dictionary of PATH priors
        multi (bool): Whether to run the simulation in multiprocessing mode
        ncpu (int): The number of CPUs to use
        debug (bool): Whether to run the simulation in debug mode

    Returns:
        pandas.DataFrame: The simulation results dataframe
    """

    # FRBs 
    if debug:
        nFRB = 100
    else:
        nFRB = len(frbs)

    # Generate the FRB dicts
    FRB_dicts = []
    maxx_box = 0.
    for index, row in frbs.iterrows():
        idict = {}
        # Localization
        idict['ra'] = row.ra
        idict['dec'] = row.dec
        idict['ltype'] = 'eellipse'
        idict['lparam'] = {'a': row.a,
                'b': row.b, 'theta': row.PA}
        # Prior
        idict['priors'] = prior_dict
        idict['index'] = index

        # Choose "local" or "fixed" grid likelihood calculations
        # (most of the time "local" will be the right answer here)
        idict['pmode'] = 'local'
        # Set the step_size method: 
        # 'relative' -- Step size is relative to the galaxy size
        # 'absolute' -- Step size is absolute in arcsec [not recommended]
        idict['step_size_mode'] = 'relative'
        # step_size should be set to 0.05 for local runs
        idict['step_size'] = 0.05

        # Box sizes
        ssize, max_box = set_anly_sizes(idict['ltype'], 
                                        idict['lparam'])
        idict['ssize'] = ssize
        idict['max_box'] = max_box
        if max_box > maxx_box:
            maxx_box = max_box

        # Save
        FRB_dicts.append(idict)

    # ####################################################
    # Build the galaxy catalog tables
    frb_coords = SkyCoord(ra=frbs.ra.values,
                          dec=frbs.dec.values, unit='deg')   
    galaxy_coords = SkyCoord(ra=catalog.ra.values,
                             dec=catalog.dec.values, unit='deg')

    # Search
    idx1, idx2, sep2d, _ = search_around_sky(
        galaxy_coords, frb_coords, maxx_box*units.arcsec)


    print("Slicing...")
    list_candidates = []
    for kk in range(nFRB):
        if (kk % 1000) == 0:
            print('kk: ', kk)
        in_idx2 = np.where(idx2 == kk)[0]
        gd_gal = idx1[in_idx2]
        close_galaxies = catalog.iloc[gd_gal][
            ['ang_size', 'mag', 'ra', 'dec', 'ID']].copy()
        # close_galaxies['separation'] = sep2d[in_idx2].to('arcsec').value
        list_candidates.append(close_galaxies)

    # Slicing done — the full catalog and coords are no longer needed
    # Delete them BEFORE creating the pool so workers don't inherit 31.5G
    del catalog
    del galaxy_coords
    del idx1, idx2, sep2d
    gc.collect()

    # Now create the pool - workers will fork from a much smaller parent
    # Run PATH
    print("PATH time")
    print("Will take a while, ~1 hr for 10,000 FRBs, depending on your computing setup and ncpu.")
    if multi:
        chunksize = max(1, nFRB // (ncpu * 8))
        args_iter = ((i, FRB_dicts[i], list_candidates[i]) for i in range(nFRB))
        all_tbls = []
        with multiprocessing.Pool(processes=ncpu) as pool:
            for ii, (idx, sv_tbl) in enumerate(
                    pool.imap_unordered(run_dict_wrapper, args_iter, chunksize=chunksize)):
                if (ii % 1000) == 0:
                    print(f'Collecting result {ii}/{nFRB}...')
                if sv_tbl is not None:
                    sv_tbl['iFRB'] = idx
                    all_tbls.append(sv_tbl)
        final_tbl = pandas.concat(all_tbls)
    else:
        all_tbls = []
        for ii in range(nFRB):
            idx, sv_tbl = run_dict_wrapper((ii, FRB_dicts[ii], list_candidates[ii]))
            if sv_tbl is not None:
                sv_tbl['iFRB'] = idx
                all_tbls.append(sv_tbl)
        final_tbl = pandas.concat(all_tbls)
    # Finish
    final_tbl.reset_index(inplace=True, drop=True)

    # Return
    return final_tbl
