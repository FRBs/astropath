import os
from pathlib import Path
import pandas
import numpy as np
import random
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units
#from zdm.chime import grids
#from zdm.loading import load_CHIME
#cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
from astropy.cosmology.realizations import Planck18 as cosmo
from astropy import units as u

from astropath.simulations import assign_frbs_to_hosts, load_galaxy_catalog

from IPython import embed

# Set the typical size of the localization region
localization = (15., 2., 12.)  # a (arcsec), b (arcsec), PA (deg)

def assign_chime_frbs_to_hosts_exponential(
    frb_tbl:pandas.DataFrame, 
    galaxy_df:pandas.DataFrame,
    outfile:str,
    localization:tuple,
    trim_catalog:units.Quantity=60*units.arcmin,
    return_df=True,
    seed = 42,
    scale:float=0.5, # half-light
    ):
    """
    Assigns CHIME FRBs to hosts based on given inputs.

    Writes the table to disk as a parquet file

    Args:
        frb_file (str): Path to the file containing FRB data.
        galaxy_catalog (pandas.DataFrame): DataFrame containing galaxy catalog data.
        outfile (str): Path to the output file where the results will be saved.
        localization (tuple): Tuple containing localization information (a, b, PA).

    """
    # Set the random seed
    random.seed(seed)
    np.random.seed(seed)

    # Prep
    nsample = len(frb_tbl['m_r'])

    ra_deg = galaxy_df.ra.values * units.deg
    dec_deg = galaxy_df.dec.values * units.deg
    ra_min, ra_max = ra_deg.min(), ra_deg.max()
    dec_min, dec_max = dec_deg.min(), dec_deg.max()

    print('Trim galaxies')
    cut_ra = (ra_deg > (ra_min + trim_catalog)) & (
        ra_deg < (ra_max - trim_catalog))
    cut_dec = (dec_deg > (dec_min + trim_catalog)) & (
        dec_deg < (dec_max - trim_catalog))

    # Cut me
    cuts = cut_dec & cut_ra
    galaxy_cut = galaxy_df[cuts]

    print('Set up galaxy magnitude lists')
    fake_coords = SkyCoord(ra=np.ones(nsample),
                           dec=frb_tbl['m_r'], unit='deg')
    fake_galaxy = SkyCoord(ra=np.ones(len(galaxy_cut)),
                           dec=galaxy_cut.mag_best.values,
                           unit='deg')

    # Prep for matching
    galaxy_flag = np.ones(len(galaxy_cut), dtype=bool)
    galaxy_flag_idx = np.arange(len(galaxy_cut))
    galaxy_idx = galaxy_cut.index.values.copy()
    frb_idx = -1*np.ones(len(fake_coords), dtype=int)

    print('Match FRBs to galaxies by magnitude')
    while(np.any(frb_idx < 0)):

        print(f"Remaining: {np.sum(frb_idx < 0)}")
        print(f"Brightest: {np.min(fake_coords.dec)}")

        # Sub me
        sub_fake_coords = fake_coords[frb_idx < 0]
        sub_frb_idx = np.where(frb_idx < 0)[0]

        sub_fake_galaxy = fake_galaxy[galaxy_flag]
        sub_galaxy_idx = galaxy_idx[galaxy_flag] # Index in the full galaxy table
        sub_galaxy_flag_idx = galaxy_flag_idx[galaxy_flag] # Index for the flagging

        # Ran out of bright ones?
        if np.max(sub_fake_coords.dec.deg) < np.min(sub_fake_galaxy.dec.deg):
            srt_galaxy = np.argsort(sub_fake_galaxy.dec.deg)
            srt_frb = np.argsort(sub_fake_coords.dec.deg)
            # Set
            frb_idx[sub_frb_idx[srt_frb]] = sub_galaxy_idx[srt_galaxy[:len(srt_frb)]]
            assert np.all(frb_idx >= 0)
            mag_bright_cut = sub_fake_galaxy.dec.deg[srt_galaxy][len(sub_fake_coords)]
            print(f'Ran out of bright ones at {mag_bright_cut}')
            #embed(header='monte_carlo.py: 153')
            break


        print(f"Min: {sub_fake_coords.dec.min()}")
        print(f"Max: {sub_fake_coords.dec.max()}")

        # Match
        idx, d2d, _ = match_coordinates_sky(
            sub_fake_coords, sub_fake_galaxy,
            nthneighbor=1)

        # Worst case
        imx = np.argmax(d2d)
        #print(f'Max: {sub_fake_coords[imx]}')
        print(f'sep = {d2d[imx]}')

        # Take a cosmo galaxy only once
        uni, uni_idx = np.unique(idx, return_index=True)

        frb_idx[sub_frb_idx[uni_idx]] = sub_galaxy_idx[uni]
        galaxy_flag[sub_galaxy_flag_idx[uni]] = False

        #if debug:
        #    imn = np.argmin(fake_coords.dec)
        #    embed(header='monte_carlo.py: 116')

    print('Generating the FRB coordinates')
    galaxy_sample = galaxy_cut.loc[frb_idx]
    galaxy_coords = SkyCoord(ra=galaxy_sample.ra.values,
                             dec=galaxy_sample.dec.values,
                             unit='deg')

    # Offset the FRB in the galaxy
    # ORIGINAL:
    # theta_max = galaxy_sample.half_light.values / scale
    # randn = np.random.normal(size=10*nsample)
    # gd = np.abs(randn) < (6.*scale)
    # randn = randn[gd][0:nsample]
    # galaxy_offset = randn * theta_max * units.arcsec
    # gal_pa = np.random.uniform(size=nsample, low=0., high=360.)
    
    # UNIFORM distribution
    # randn = np.random.uniform(low=0., high=10., size=nsample)
    # galaxy_offset = randn * galaxy_sample.half_light.values * units.arcsec
    # gal_pa = np.random.uniform(size=nsample, low=0., high=360.)
    
    print('EXPONENTIAL distribution')
    randn = np.random.exponential(scale=scale, size=10*nsample)
    gd = np.abs(randn) < (6.)
    randn = randn[gd][0:nsample]
    galaxy_offset = randn * galaxy_sample.half_light.values * units.arcsec
    gal_pa = np.random.uniform(size=nsample, low=0., high=360.)

    embed(header='test_sim_assign 155')
    print("Offsetting FRBs in galaxy...")
    frb_coord = [coord.directional_offset_by(
        gal_pa[kk]*units.deg, galaxy_offset[kk]) 
                 for kk, coord in enumerate(galaxy_coords)]

    true_frb_coord = frb_coord.copy()

    # Offset by Localization
    randn = np.random.normal(size=10*nsample)
    gd = np.abs(randn) < 3. # limit to 3 sigma
    randn = randn[gd]

    # TODO -- Make sure this is right
    a_off = randn[0:nsample] * localization[0] * units.arcsec
    b_off = randn[nsample:2*nsample] * localization[1] * units.arcsec
    #pa = np.arctan2(decoff, raoff) * 180./np.pi - 90.
    local_offset = np.sqrt(a_off**2 + b_off**2) 

    print("Offsetting FRB by localization...")
    frb_coord = [coord.directional_offset_by(
        localization[2]*units.deg, a_off[kk])
                 for kk, coord in enumerate(frb_coord)]
    frb_coord = [coord.directional_offset_by(
        (localization[2]+90)*units.deg, b_off[kk])
                 for kk, coord in enumerate(frb_coord)]

    # Write to disk
    df = pandas.DataFrame()
    df['ra'] = [coord.ra.deg for coord in frb_coord]
    df['dec'] = [coord.dec.deg for coord in frb_coord]
    df['true_ra'] = [coord.ra.deg for coord in true_frb_coord]
    df['true_dec'] = [coord.dec.deg for coord in true_frb_coord]
    df['gal_ID'] = galaxy_sample.ID.values
    df['gal_off'] = galaxy_offset.value # arcsec
    df['mag'] = galaxy_sample.mag_best.values
    df['half_light'] = galaxy_sample.half_light.values
    df['loc_off'] = local_offset.value # arcsec

    # Add ID
    df['FRB_ID'] = np.arange(len(frb_tbl))

    # Add localization
    df['a'] = localization[0]
    df['b'] = localization[1]
    df['PA'] = localization[2]

    df.to_parquet(outfile, index=False)
    print(f"Wrote: {outfile}")
    
    if return_df:
        return df


def test_orig_assign(seed:int=42, NFRB:int=100):

    
    # Big catalog
    combined_file = os.path.join(os.getenv('FRB_APATH'), 'combined_HSC_DECaLs_HECATE_galaxies_hecatecut.parquet')
    catalog = pandas.read_parquet(combined_file)

    # Load up the generated FRBs
    generated_frbs_fn = os.path.join('generated_frbs_test_{}_{}.parquet'.format(int(seed), int(NFRB)))
    frbs = pandas.read_parquet(generated_frbs_fn)
    
    # Output file
    output_fn = 'generated_hosts_DECaL_DECaLhost_hecatecut_{}_{}.parquet'.format(int(seed), int(NFRB))

    # Assign
    hosts = assign_chime_frbs_to_hosts_exponential(
        frbs, catalog, output_fn, localization)

def test_astropath_assign(seed:int=42, NFRB:int=100):

    # Load FRBs
    output_fn = 'generated_frbs_test_{}_{}.parquet'.format(int(seed), int(NFRB))
    df_frbs = pandas.read_parquet(output_fn)

    # Load real catalog 
    galaxies = load_galaxy_catalog()

    # Assign
    hosts = assign_frbs_to_hosts(
        df_frbs, galaxies, localization, seed=seed,
        trim_catalog=60*units.arcmin)

    # Load original
    output_fn = 'generated_hosts_DECaL_DECaLhost_hecatecut_{}_{}.parquet'.format(int(seed), int(NFRB))
    hosts_orig = pandas.read_parquet(output_fn)

    #embed(header='test_astropath_assign 239')
    # Test
    assert np.allclose(hosts['gal_ID'].values, hosts_orig['gal_ID'].values)
    assert np.allclose(hosts['ra'].values, hosts_orig['ra'].values)
    assert np.allclose(hosts['dec'].values, hosts_orig['dec'].values)
    assert np.allclose(hosts['true_ra'].values, hosts_orig['true_ra'].values)
    assert np.allclose(hosts['true_dec'].values, hosts_orig['true_dec'].values)
    assert np.allclose(hosts['gal_off'].values, hosts_orig['gal_off'].values)
    assert np.allclose(hosts['mag'].values, hosts_orig['mag'].values)
    assert np.allclose(hosts['half_light'].values, hosts_orig['half_light'].values)

    print("Tests passed!!")

# Command line
if __name__ == '__main__':
    import sys
    flg = int(sys.argv[1])

    if flg == 0:
        test_orig_assign()
    elif flg == 1:
        test_astropath_assign()
    else:
        raise ValueError(f"Invalid flag: {flg}")
    #test_astropath_assign()