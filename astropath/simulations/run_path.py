
import numpy as np
import pandas
import multiprocessing

from astropy.coordinates import SkyCoord
from astropy import units
from astropy.coordinates import search_around_sky
from astropy.table import Table

from astropath.run import run_on_dict, set_anly_sizes

from IPython import embed

    

def run_dict_wrapper(idx:int, idict:dict,
                     catalog:pandas.DataFrame):
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
    # Convert to astropy table
    catalog = Table.from_pandas(catalog)
    # Run me
    candidates, P_Ux, Path, mag_key, cut_catalog, stars = \
        run_on_dict(idict, catalog=catalog, mag_key='mag')

    # Return
    sv_tbl = candidates.sort_values(
        'P_Ox', ascending=False)
    sv_tbl['gal_ID'] = sv_tbl.index.values

    # Return
    return sv_tbl, idx, Path
    

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
        # 
        idict = {}
        # Localization
        idict['ra'] = row.ra
        idict['dec'] = row.dec
        idict['ltype'] = 'eellipse'
        idict['lparam'] = {'a': row.a,
                'b': row.b, 'theta': row.PA}
        # Prior
        idict['priors'] = prior_dict

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
        # Grab em
        in_idx2 = np.where(idx2 == kk)[0]
        gd_gal = idx1[in_idx2]
        close_galaxies = catalog.iloc[gd_gal][
            ['ang_size', 'mag', 'ra', 'dec', 'ID']]
        # Extras
        close_galaxies['separation'] = sep2d[in_idx2].to('arcsec').value
        close_galaxies['coords'] = galaxy_coords[gd_gal]
        #
        list_candidates.append(close_galaxies)

    # RUn it!
    print("PATH time")
    if multi:
        pool = multiprocessing.Pool(processes=ncpu)
        results = [pool.apply_async(run_dict_wrapper,
            args=(idx_FRB, FRB_dicts[idx_FRB], 
                  list_candidates[idx_FRB]))
            for idx_FRB in range(nFRB)]
        # Run
        output = [p.get() for p in results]
        idx = [item[1] for item in output]
        # Unpack
        all_tbls = np.array([item[0] for item in output], dtype=object)
        all_tbls = all_tbls[idx]
        # Expunge the None's
        gd_tbl = np.array([False if item is None else True for item in all_tbls])
        gd_idx = np.arange(all_tbls.size)[gd_tbl]
        all_tbls = all_tbls[gd_tbl]
        for kk in range(all_tbls.size):
            all_tbls[kk]['iFRB'] = gd_idx[kk]

    else:
        idx_FRB = 0
        results = run_dict_wrapper(
            idx_FRB, FRB_dicts[idx_FRB],
            list_candidates[idx_FRB])

    # Finish
    final_tbl = pandas.concat(all_tbls)
    final_tbl.reset_index(inplace=True, drop=True)

    # Return
    return final_tbl