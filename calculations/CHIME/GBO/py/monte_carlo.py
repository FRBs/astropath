""" Monte Carlo calculations for CHIME/FRB with GBO """
import numpy as np

from scipy import interpolate

import pandas

from matplotlib import pyplot as plt
import seaborn as sns

from astropy import units

from astropy.coordinates import SkyCoord, match_coordinates_sky

from astropath import cosmos
from astropath import montecarlo

from IPython import embed

# Localization error assumed for GBO
gbo_radec_sigma=(3.,15.) # arcsec, ra,dec

def generate_frbs(outfile:str,
    chime_mr_tbl:str='CHIME_mr_5Jyms_150.parquet', 
                  m_r_min:float=15.,
                  m_r_max:float=26.,
                  radec_sigma:tuple=None,
                  scale:float=2., # half-light
                  nsample=10000,
                  debug:bool=False,
                  plots:bool=False):
    """ Genereate a set of FRBs associated to COSMOS galaxies
    With random placement in the galaxy
    and a random location error

    A pandas table is generated and saved to disk

    Args:
        outfile (str): 
            Filename for the output table (pandas)
        chime_mr_tbl (str, optional): Defaults to 'CHIME_mr_5Jyms_150.parquet'.
            Table of estimated CHIME/FRB m_r values for FRBs
        m_r_min (float, optional): Defaults to 15..
            Minimum m_r value for the CHIME/FRB distribution
        m_r_max (float, optional): Defaults to 26..
            Maximum m_r value for the CHIME/FRB distribution
        radec_sigma (tuple, optional): Defaults to (3.,15.).
            Uncertainty in the RA,DEC of the FRB in arcsec
            CHIME/GBO will be precise in RA
        scale (float, optional): Defaults to 2..
            Scale factor for the galaxy half-light radius
        debug (bool, optional): Defaults to False.
            Debugging mode
        plots (bool, optional): Defaults to False.
            Generate plots
    """
    if radec_sigma is None:
        radec_sigma = gbo_radec_sigma

    # Load up m_r distribution
    chime_mr = pandas.read_parquet(chime_mr_tbl)
    m_r_cuts = (chime_mr.m_r < m_r_max) & (
            chime_mr.m_r > m_r_min)
    chime_mr = chime_mr[m_r_cuts]
   
    if debug:
        chime_mr=chime_mr[0:100000]
    n_chime = len(chime_mr)

    # Generate a CDF
    srt = np.argsort(chime_mr.m_r.values)

    f_cdf = interpolate.interp1d(
        np.arange(n_chime)/(n_chime-1),
        chime_mr.m_r.values[srt])

    # Random sample
    rand = np.random.uniform(size=nsample)
    rand_mr = f_cdf(rand)

    if plots and False:
        plt.clf()
        ax = plt.gca()

        # Full
        cdf = np.linspace(0., 1., 10000)
        mr = f_cdf(cdf)
        ax.plot(mr, cdf, 'k-', label='Full CDF')

        # Realization
        srt2 = np.argsort(rand_mr)
        ax.step(rand_mr[srt2], 
                np.arange(nsample)/(nsample-1),
                color='r')

        ax.legend()
        plt.show()

    # Load COSMOS
    cosmos_df = cosmos.load_galaxies()
    ra_deg = cosmos_df.ra.values * units.deg
    dec_deg = cosmos_df.dec.values * units.deg
    ra_min, ra_max = ra_deg.min(), ra_deg.max()
    dec_min, dec_max = dec_deg.min(), dec_deg.max()

    cut_ra = (ra_deg > (ra_min + 1*units.arcmin)) & (
        ra_deg < (ra_max - 1*units.arcmin))
    cut_dec = (dec_deg > (dec_min + 1*units.arcmin)) & (
        dec_deg < (dec_max - 1*units.arcmin))

    # Cut me
    cuts = cut_dec & cut_ra
    cosmos_cut = cosmos_df[cuts]

    # Choose the galaxies
    fake_coords = SkyCoord(ra=np.ones(nsample),
                           dec=rand_mr, unit='deg')
    fake_cosmos = SkyCoord(ra=np.ones(len(cosmos_cut)),
                           dec=cosmos_cut.mag_best.values,
                           unit='deg')

    cosmos_flag = np.ones(len(cosmos_cut), dtype=bool)
    cosmos_flag_idx = np.arange(len(cosmos_cut))
    cosmos_idx = cosmos_cut.index.values.copy()
    frb_idx = -1*np.ones(len(fake_coords), dtype=int)

    while(np.any(frb_idx < 0)):

        print(f"Remaining: {np.sum(frb_idx < 0)}")
        # Sub me
        sub_fake_coords = fake_coords[frb_idx < 0]
        sub_frb_idx = np.where(frb_idx < 0)[0]

        sub_fake_cosmos = fake_cosmos[cosmos_flag]
        sub_cosmos_idx = cosmos_idx[cosmos_flag] # Index in the full cosmos table
        sub_cosmos_flag_idx = cosmos_flag_idx[cosmos_flag] # Index for the flagging

        # Ran out of bright ones?
        if np.max(sub_fake_coords.dec.deg) < np.min(sub_fake_cosmos.dec.deg):
            srt_cosmos = np.argsort(sub_fake_cosmos.dec.deg)
            srt_frb = np.argsort(sub_fake_coords.dec.deg)
            # Set
            frb_idx[sub_frb_idx[srt_frb]] = sub_cosmos_idx[srt_cosmos[:len(srt_frb)]]
            assert np.all(frb_idx >= 0)
            break


        print(f"Min: {sub_fake_coords.dec.min()}")
        print(f"Max: {sub_fake_coords.dec.max()}")

        # Match
        idx, d2d, _ = match_coordinates_sky(
            sub_fake_coords, sub_fake_cosmos,
            nthneighbor=1)

        # Worst case
        imx = np.argmax(d2d)
        #print(f'Max: {sub_fake_coords[imx]}')
        print(f'sep = {d2d[imx]}')

        # Take a cosmo galaxy only once
        uni, uni_idx = np.unique(idx, return_index=True)

        frb_idx[sub_frb_idx[uni_idx]] = sub_cosmos_idx[uni]
        cosmos_flag[sub_cosmos_flag_idx[uni]] = False

        #if debug:
        #    imn = np.argmin(fake_coords.dec)
        #    embed(header='monte_carlo.py: 116')

    # Generating the FRB coordinates
    cosmos_sample = cosmos_cut.loc[frb_idx]
    galaxy_coords = SkyCoord(ra=cosmos_sample.ra.values,
                             dec=cosmos_sample.dec.values,
                             unit='deg')

    # Offset the FRB in the galaxy
    theta_max = cosmos_sample.half_light.values / scale
    randn = np.random.normal(size=10*nsample)
    gd = np.abs(randn) < (6.*scale)
    randn = randn[gd][0:nsample]
    galaxy_offset = randn * theta_max * units.arcsec
    gal_pa = np.random.uniform(size=nsample, low=0., high=360.)

    print("Offsetting FRB in galaxy...")
    frb_coord = [coord.directional_offset_by(
        gal_pa[kk]*units.deg, galaxy_offset[kk]) 
                 for kk, coord in enumerate(galaxy_coords)]

    # Offset by Localization
    randn = np.random.normal(size=10*nsample)
    gd = np.abs(randn) < 3. # limit to 3 sigma
    randn = randn[gd]

    raoff = randn[0:nsample] * radec_sigma[0]
    decoff = randn[nsample:2*nsample] * radec_sigma[1]
    pa = np.arctan2(decoff, raoff) * 180./np.pi - 90.
    local_offset = np.sqrt(raoff**2 + decoff**2) * units.arcsec

    if plots:
        sns.histplot(x=pa)
        plt.show()

    print("Offsetting FRB by localization...")
    frb_coord = [coord.directional_offset_by(
        pa[kk]*units.deg, local_offset[kk]) 
                 for kk, coord in enumerate(frb_coord)]

    # Write to disk
    df = pandas.DataFrame()
    df['ra'] = [coord.ra.deg for coord in frb_coord]
    df['dec'] = [coord.dec.deg for coord in frb_coord]
    df['gal_ID'] = cosmos_sample.index.values
    df['gal_off'] = galaxy_offset.value # arcsec
    df['loc_off'] = local_offset.value # arcsec

    df.to_csv(outfile, index=False)
    print(f"Wrote: {outfile}")

def run_mc(outfile:str, debug:bool=False):
    # Load up
    frb_tbl = pandas.read_csv('frb_monte_carlo.csv')
    cosmos_df = cosmos.load_galaxies()


    if debug:
        frb_tbl = frb_tbl.iloc[0:30]

    # Define items
    path_prior = dict(U=0.1, # P_U
                      S=2.0) # scale factor
    frb_ee = dict(a=gbo_radec_sigma[1], 
                  b=gbo_radec_sigma[0],
                  theta=0.)

    # Run
    #multi = False if debug else True
    multi = True
    montecarlo.run_em(frb_tbl, cosmos_df, path_prior,
                      outfile, frb_ee, box_hwidth=60.,
                      step_size=0.2, 
                      mag_lim=23.2, # PS2
                      debug=debug, multi=multi)

# Command line execution
if __name__ == '__main__':

    # Generate FRBs
    #generate_frbs('frb_monte_carlo.csv',
    #    debug=True, plots=False, nsample=10000)

    # Monte Carlo
    run_mc('first_try.csv', debug=True)