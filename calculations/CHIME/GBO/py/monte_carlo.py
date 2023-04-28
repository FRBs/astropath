""" Monte Carlo calculations for CHIME/FRB with GBO """
import numpy as np

from scipy import interpolate

import pandas

from matplotlib import pyplot as plt
import seaborn as sns

from astropy import units

from astropath import cosmos
from astropy.coordinates import SkyCoord, match_coordinates_sky

from IPython import embed

def generate_frbs(chime_mr_tbl:str='CHIME_mr_5Jyms_150.parquet', 
                  m_r_min:float=15.,
                  m_r_max:float=26.,
                  nsample=10000,
                  debug:bool=False,
                  plots:bool=False):

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

    if plots:
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
    cosmos_idx = cosmos_cut.index.values.copy()
    frb_idx = -1*np.ones(len(fake_coords), dtype=int)

    while(np.any(frb_idx < 0)):

        print(f"Remaining: {np.sum(frb_idx < 0)}")
        # Sub me
        sub_fake_coords = fake_coords[frb_idx < 0]
        sub_frb_idx = np.where(frb_idx < 0)[0]
        sub_fake_cosmos = fake_cosmos[cosmos_flag]
        sub_cosmos_idx = cosmos_idx[cosmos_flag]

        print(f"Min: {sub_fake_coords.dec.min()}")
        print(f"Max: {sub_fake_coords.dec.max()}")

        # Match
        idx, d2d, _ = match_coordinates_sky(
            sub_fake_coords, sub_fake_cosmos,
            nthneighbor=1)

        imx = np.argmax(d2d)
        #print(f'Max: {sub_fake_coords[imx]}')
        print(f'sep = {d2d[imx]}')

        # Specify them
        uni, uni_idx = np.unique(idx, return_index=True)

        if debug:
            embed(header='monte_carlo.py: 116')
        cosmos_flag[sub_cosmos_idx[uni_idx]] = False
        frb_idx[sub_frb_idx[uni_idx]] = sub_cosmos_idx[idx[uni_idx]]

    embed(header='monte_carlo.py: 113')


# Command line execution
if __name__ == '__main__':
    generate_frbs(debug=True, plots=False)