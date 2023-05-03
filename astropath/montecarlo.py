""" Methods to perform Monte Carlo simlulations """

import sys, os
import numpy as np
import pandas
import multiprocessing

from astropy import units
from astropy.coordinates import SkyCoord
from astropy.coordinates import search_around_sky

from frb.frb import FRB

from astropath import path 
from astropath import cosmos

from IPython import embed

def path_calc(idx:int, frb_coord:SkyCoord, frb_ee:dict, 
              close_galaxies:pandas.DataFrame, 
              prior:dict, box_hwidth:float,
              step_size:float):
    """
    Perform the main PATH calculation

    Args:
        idx (int): Simple index for referencing when multiprocessing
        frb_coord (SkyCoord): FRB sky coordinate
        frb_ee (dict): FRB error ellipse parameters
        close_galaxies (pandas.DataFrame): 
            Table of galaxies near the FRB
        prior (dict): _description_
        box_hwidth (float): 
            Half-width of the box to analyze galaxies (arcsec)
        step_size (float):
            Step size for the grid analysis (arcsec)
            0.2" is ok for COSMOS 

    Returns:
        tuple:  pandas.DataFrame, int, path.PATH
    """


    # Empty?
    if len(close_galaxies) == 0:
        return None, idx, None

    # Turn into a cndidate table
    Path = path.PATH()

   # Set up localization
    Path.init_localization('eellipse', 
                           center_coord=frb_coord, 
                           eellipse=frb_ee.copy())
    # Candidates
    Path.init_candidates(close_galaxies['ra'],
                         close_galaxies['dec'],
                         close_galaxies['ang_size'],
                         mag=close_galaxies['mag'])

    # Candidate prior
    Path.init_cand_prior('inverse', P_U=prior['U'])

    # Offset prior
    Path.init_theta_prior('exp', 6., prior['S'])

    # Priors
    p_O = Path.calc_priors()

    # Posterior
    #embed(header='68 of montecarlo.py')
    P_Ox, P_Ux = Path.calc_posteriors(
        'fixed', box_hwidth=box_hwidth, 
        step_size=step_size, 
        max_radius=box_hwidth)

    # Save
    sv_tbl = Path.candidates.sort_values(
        'P_Ox', ascending=False)#.drop(columns=['coords'])

    # Return
    return sv_tbl, idx, Path


def run_em(mc_frbs:pandas.DataFrame, mc_galaxies:pandas.DataFrame, 
           path_prior:dict, outfile:str, frb_ee:dict, 
           galaxy_coords=None, 
           box_hwidth:float=20., step_size:float=0.1,
           mag_lim:float=None,
           ncpu:int=15, multi:bool=True, mc_defs:dict=None,
           debug:bool=False):
    # Defs
    if mc_defs is None:
        mc_defs = cosmos.cosmos_defs()

    # Load FRBs
    print("Loading...")
    nFRB = len(mc_frbs)

    # frb_ee
    if isinstance(frb_ee, dict):
        frb_ee_list = [frb_ee]*nFRB
    elif isinstance(frb_ee, str):
        frb_ee_list = []
        if frb_ee == 'loc_sig':
            for kk in range(nFRB):
                frb_ee_list.append({'a': mc_frbs.iloc[kk].loc_sig,
                                    'b': mc_frbs.iloc[kk].loc_sig,
                                    'theta': 0.})
    else:
        raise IOError("Bad frb_ee")

    # Cut
    mc_galaxies = mc_galaxies[np.isfinite(mc_galaxies.a_image) &
                              np.isfinite(mc_galaxies.mag_best)]

    # Extra bits
    mc_galaxies['ang_size'] = mc_galaxies.a_image * mc_defs['plate_scale']
    #mc_galaxies[mc_defs['filter']] = mc_galaxies.mag_best
    mc_galaxies['mag'] = mc_galaxies.mag_best

    # Cut on mag?
    if mag_lim is not None:
        keep = mc_galaxies.mag_best < mag_lim
        mc_galaxies = mc_galaxies[keep]

    # Coordinates (this is slow)
    print("Building coords...")
    frb_coords = SkyCoord(ra=mc_frbs.ra, dec=mc_frbs.dec, unit='deg')
    if galaxy_coords is None:
        galaxy_coords = SkyCoord(ra=mc_galaxies.ra, dec=mc_galaxies.dec, unit='deg')

    # Find the candidates
    print("Finding candidates...")
    idx1, idx2, sep2d, _ = search_around_sky(
        galaxy_coords, frb_coords, box_hwidth*units.arcsec)

    # Slice the tables
    print("Slicing...")
    list_candidates = []

    #nFRB = 2000
    for kk in range(nFRB):
        if (kk % 1000) == 0:
            print('kk: ', kk)
        # Grab em
        in_idx2 = np.where(idx2 == kk)[0]
        gd_gal = idx1[in_idx2]
        close_galaxies = mc_galaxies.iloc[gd_gal][
            ['ang_size', 'mag', 'ra', 'dec']]
        # Extras
        close_galaxies['separation'] = sep2d[in_idx2].to('arcsec').value
        close_galaxies['coords'] = galaxy_coords[gd_gal]
        #
        list_candidates.append(close_galaxies)

    # Loop me!
    print('Starting to loop!')
    # Build
    if multi:
        pool = multiprocessing.Pool(processes=ncpu)
        results = [pool.apply_async(path_calc,
            args=(idx_FRB, frb_coords[idx_FRB], frb_ee_list[idx_FRB],
                    list_candidates[idx_FRB], path_prior, 
                    box_hwidth, step_size))
            for idx_FRB in range(nFRB)]
        output = [p.get() for p in results]
        idx = [item[1] for item in output]
        all_tbls = np.array([item[0] for item in output], dtype=object)
        all_tbls = all_tbls[idx]
        # Expunge the None's
        gd_tbl = np.array([False if item is None else True for item in all_tbls])
        gd_idx = np.arange(all_tbls.size)[gd_tbl]
        all_tbls = all_tbls[gd_tbl]
        for kk in range(all_tbls.size):
            all_tbls[kk]['iFRB'] = gd_idx[kk]
        #embed(header='montecarlo: 177')

    else:
        idx_FRB = 0
        results = path_calc(idx_FRB, frb_coords[idx_FRB], frb_ee_list[idx_FRB],
                    list_candidates[idx_FRB], path_prior, 
                    box_hwidth, step_size)
        embed(header='183 of montecarlo.py')
        '''
        all_tbls = []
        for idx_FRB in range(nFRB):
            sv_tbl = path_calc(frb_coords[idx_FRB], mc_galaxies, galaxy_coords,
                      FRB, frbA, mc_defs, prior)

            sv_tbl['iFRB'] = idx_FRB
            all_tbls.append(sv_tbl)
        '''
    if debug:
        embed(header='montecarlo: 195')

    # Finish
    final_tbl = pandas.concat(all_tbls)
    final_tbl.to_csv(outfile)
    print("Wrote: {}".format(outfile))

