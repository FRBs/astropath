from importlib.resources import files

import pandas
import numpy as np

import scipy.stats as stats
from frb.galaxies import hosts as hosts_mod
from frb.frb_surveys import chime
from scipy.interpolate import interp1d
from frb.defs import frb_cosmo
from astropy.cosmology.realizations import Planck18
from astropy import units

# JXP imports
from frb.dm import prob_dmz

# astropath
from astropath.simulations import generate_frbs, SURVEY_GRIDS
from astropath.simulations.generate_frbs import load_chime_cat1_DMeg

from IPython import embed

import random

def orig_generate():
    # Set the random seed
    seed = 42
    random.seed(seed)
    np.random.seed(seed)


    NFRB = 100

    # Load p(z|DM) grid
    chime_grid = prob_dmz.grab_repo_grid('CHIME_pzdm.npz')
    zvals = chime_grid['z']
    dmvals = chime_grid['DM']
    all_rates = chime_grid['pzdm']

    # Load CHIME Catalog 1
    DMex = load_chime_cat1_DMeg()

    # Cumulative
    cum_all = np.cumsum(all_rates, axis=0)
    norm = np.outer(np.ones(zvals.size), cum_all[-1,:])
    cum_all /= norm
    cum_all[0,:] = 0.

    # Interpolators
    print("Building interpolators")
    fs = [interp1d(cum_all[:,ii], zvals) for ii in range(dmvals.size)]

    kernel = stats.gaussian_kde(DMex)#, bw_method=0.6)
    dms = np.linspace(0., 3000, 500)
    DMex_kde = kernel(dms)

    cum_DMex = np.cumsum(DMex_kde)
    cum_DMex[0] = 0.
    fh = interp1d(cum_DMex/cum_DMex[-1], dms)

    # Random numbers
    rand = np.random.uniform(size=NFRB)
    rand_DMex = fh(rand) # np.random.choice(df_dr1['DMex'], size=NFRB) # 

    # Grab redshifts
    zs = []
    rand = np.random.uniform(size=NFRB)
    for kk,DMc in enumerate(rand_DMex):
        imin = np.argmin(np.abs(dmvals-DMc))
        z = fs[imin](rand[kk])
        zs.append(float(z))
    zs = np.array(zs)

    # Original Mr PDF
    orig_Mr = False
    if orig_Mr:
        # Load the previous distribution
        Mr, density = hosts_mod.load_Mr_pdf()

        cum_Mr = np.cumsum(density)
        cum_Mr[0] = 0.
        fMr = interp1d(cum_Mr/cum_Mr[-1], Mr)
    else:
        # Load up Lz values
        host_file = files('astropath.data') / 'frb_surveys' / 'Lz_host_data.csv'
        df= pandas.read_csv(host_file)

        # Scale mrs with z's to find distribution of Mrs
        mrs = np.array(df['r-band'])
        zs_mrs = np.array(df['redshift'])

        # Get luminosity distance
        ds = Planck18.luminosity_distance(zs_mrs).to(units.parsec).value

        # Calculate absolute magnitudes
        Mrs = mrs - 5. * np.log10(ds) + 5

        # Calculate KDE of absolute magnitudes
        kernel = stats.gaussian_kde(Mrs)#, bw_method=0.6)
        mags = np.linspace(-25., -15., 500)
        Mr_kde = kernel(mags)

        # Cum sum
        cum_Mr = np.cumsum(Mr_kde)
        cum_Mr[0] = 0.
        fMr = interp1d(cum_Mr/cum_Mr[-1], mags)

    # Go forth

    rand = np.random.uniform(size=NFRB)
    rand_Mr = fMr(rand)

    dist_mod = frb_cosmo.distmod(zs).value
    host_m_r = dist_mod + rand_Mr

    # Grab P(z,DM) grid to plot
    grid = all_rates
    z = zvals
    DM_EG = dmvals
    pzDM = grid.flatten()
    cum_sum = np.cumsum(pzDM)
    # Normalize
    cum_sum = cum_sum/cum_sum[-1]

    '''
    # Interpolate
    DM, Z = np.meshgrid(DM_EG, z)  # 2D grid for interpolation
    interp = CloughTocher2DInterpolator(list(zip(Z.ravel(), DM.ravel())), grid.ravel(), fill_value=np.nan)
    # Make new regular grid of z, dm coordinates
    z_new_list = np.linspace(0, 2.5, 500)
    dm_new_list = np.linspace(0, 4000, 500)
    z_new, dm_new = np.meshgrid(z_new_list, dm_new_list)
    Pz_dm_interp = interp(z_new, dm_new)
    Pz_dm_interp[Pz_dm_interp < 0.] = 0.
    '''

    # Build FRB table
    df_frbs = pandas.DataFrame()
    df_frbs['DMeg'] = rand_DMex
    df_frbs['z'] = zs
    df_frbs['M_r'] = rand_Mr
    df_frbs['m_r'] = host_m_r
    output_fn = 'generated_frbs_test_{}_{}.parquet'.format(int(seed), int(NFRB))
    df_frbs.to_parquet(output_fn)
    print(f"Generated {NFRB} FRBs and saved to {output_fn}")


def astropath_generate():
    # astropath version
    NFRB = 100
    seed = 42
    #random.seed(seed)
    #np.random.seed(seed)

    # Load CHIME Catalog 1
    DMeg = load_chime_cat1_DMeg()

    #df_dr1 = pandas.read_csv('./chimefrbcat1.csv')
    ## Make cut based on bonsai S/N
    #cut_snr = 12.
    #snr_cut = df_dr1['bonsai_snr'] > cut_snr
    #df_dr1 = df_dr1[snr_cut].copy()
    #DMex = np.nanmean([df_dr1['dm_exc_ne2001'].values,df_dr1['dm_exc_ymw16'].values], axis=0) #- 100. # Minus MW halo component

    # Generate
    df_chime = generate_frbs(NFRB, 'CHIME', seed=seed, 
        dm_catalog=DMeg, dm_range=(0., 3000.))

    # Load 
    output_fn = 'generated_frbs_test_{}_{}.parquet'.format(int(seed), int(NFRB))
    df_orig = pandas.read_parquet(output_fn)
    
    # Test
    assert np.allclose(df_chime['DMeg'].values, df_orig['DMeg'].values)
    assert np.allclose(df_chime['z'].values, df_orig['z'].values)
    assert np.allclose(df_chime['M_r'].values, df_orig['M_r'].values)
    assert np.allclose(df_chime['m_r'].values, df_orig['m_r'].values)
    print("Tests passed!!")


# Command line
if __name__ == '__main__':
    #orig_generate()
    astropath_generate()