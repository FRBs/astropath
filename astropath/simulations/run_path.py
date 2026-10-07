
import numpy as np
import pandas
import multiprocessing
assert multiprocessing.get_start_method() == "fork", "USR_raw_prior_Oi will not reach worker processes under this start method"

from astropy.coordinates import SkyCoord
from astropy import units
from astropy.coordinates import search_around_sky
from astropy.table import Table
from scipy.spatial import cKDTree
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
           frb_chunk:int=5000,
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

    print("Generate the FRB dicts")
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
    print("Galaxy catalog cross-match")
    # Per-FRB search radius, NOT the global maximum over all FRBs.
    #
    # run_on_dict() cuts the catalog at idict['ssize']*60, which is exactly this
    # FRB's own max_box (astropath/run.py, "Cut down the catalog based on ssize").
    # So every galaxy the old global-radius search returned beyond an FRB's own
    # max_box was thrown away before PATH ever saw it. Using max(max_box) for all
    # FRBs inflated the pair list by (max(max_box) / max_box_i)^2 while changing
    # no result: for the basecat2 localizations (p50 25", max 152") that is ~16x,
    # i.e. ~370M pairs / 33 GB at 100k FRBs against ~23M / 2 GB.
    frb_radii = np.array([d['max_box'] for d in FRB_dicts], dtype=float)  # arcsec
    print(f"  max_box: median {np.median(frb_radii[:nFRB]):.0f}\", "
        f"max {frb_radii[:nFRB].max():.0f}\"  "
        f"(the previous code used the max for every FRB)")
    
    print("Slicing...")
    cols = ['ang_size', 'mag', 'ra', 'dec', 'ID']
    col_arrays = {c: catalog[c].to_numpy() for c in cols}
    cat_index = catalog.index.to_numpy()
    
    def _unit_vec(ra_deg, dec_deg):
      """(N,3) unit vectors on the sphere.
    
      Chord length is monotonic in angular separation, so a chord cut is an
      exact angular cut. This also avoids building SkyCoord objects for the
      whole catalog, which was itself a large allocation.
      """
      ra = np.radians(np.asarray(ra_deg, dtype=float))
      dec = np.radians(np.asarray(dec_deg, dtype=float))
      cd = np.cos(dec)
      return np.column_stack([cd * np.cos(ra), cd * np.sin(ra), np.sin(dec)])
    
    tree = cKDTree(_unit_vec(catalog.ra.values, catalog.dec.values))
    
    # Chunk the FRBs so the pair list never exists for all of them at once.
    # Peak memory is set by frb_chunk, not by nFRB.
    frb_ra = frbs.ra.values
    frb_dec = frbs.dec.values
    list_candidates = []
    n_pairs = 0
    for lo in range(0, nFRB, frb_chunk):
      hi = min(lo + frb_chunk, nFRB)
      vec = _unit_vec(frb_ra[lo:hi], frb_dec[lo:hi])
      chord = 2.0 * np.sin(np.radians(frb_radii[lo:hi] / 3600.) / 2.)
      for hits in tree.query_ball_point(vec, chord):
          gd_gal = np.sort(np.asarray(hits, dtype=np.int64))
          n_pairs += len(gd_gal)
          list_candidates.append(pandas.DataFrame(
              {c: col_arrays[c][gd_gal] for c in cols},
              index=cat_index[gd_gal],
          ))
      del vec, chord
      print(f"  {hi}/{nFRB} FRBs sliced, {n_pairs:,d} candidates so far")
    
    del col_arrays, cat_index, tree
    
    # Slicing done -- the full catalog is no longer needed.
    # Delete it BEFORE creating the pool so workers do not inherit it.
    del catalog
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
