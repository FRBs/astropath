""" Figures related to CHIME/GBO calculations """
# imports
import numpy as np

import pandas
import argparse

import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Ellipse

from astropy.coordinates import SkyCoord

from frb.figures import utils as ffutils
from astropath import cosmos

import analysis

from IPython import embed

def fig_roc(path_file:str, frb_file:str,
                     outfile:str): 
    # parse
    frbs = analysis.parse_PATH(path_file, frb_file)

    # Add success
    success = frbs.PATH_ID == frbs.gal_ID
    frbs['success'] = success

    # 

    def gen_cdf(tbl):
        
        # Sort on mag
        mag_srt = np.argsort(tbl.mag.values)

        # CDF
        cdf = np.cumsum(
            tbl.success.values[mag_srt]
            )/np.cumsum(np.ones(len(tbl)))

        cut_90 = np.argmin(np.abs(cdf - 0.9))
        cut_99 = np.argmin(np.abs(cdf - 0.99))
        print(f'90% completeness at {tbl.mag.values[mag_srt][cut_90]:g}')
        print(f'99% completeness at {tbl.mag.values[mag_srt][cut_99]:g}')

        return mag_srt, cdf 

    # Figure
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111)

    # All
    mag_srt, cdf_all = gen_cdf(frbs)
    ax.plot(frbs.mag.values[mag_srt], cdf_all, 'g-', lw=2, label='All')

    # High P_Ox
    POx_90 = frbs.P_Ox > 0.9
    mag_90, cdf_90 = gen_cdf(frbs[POx_90])
    #embed(header='54 of figs')
    ax.plot(frbs[POx_90].mag.values[mag_90], cdf_90, 'b-', 
            lw=2, label='P(O|x) > 90%')

    # Label
    ax.set_xlabel('Magnitude')
    ax.set_ylabel('TP Fraction')

    ax.legend()

    # Finish
    ffutils.set_fontsize(ax, 16.)

    pad = 0.2
    plt.tight_layout(pad=0.2,h_pad=pad,w_pad=pad)

    plt.savefig(outfile, dpi=300)#, bbox_inches='tight')
    plt.close()
    print(f"Wrote: {outfile}")

    # Final stat
    gd_mag = frbs.mag < 20.5
    smpl = gd_mag & POx_90
    print(f'We might target {np.sum(smpl)}/{len(frbs)}')

def fig_scatter_Pm(path_file:str, frb_file:str,
                     outfile:str): 
    # parse
    frbs = analysis.parse_PATH(path_file, frb_file)

    ax = sns.scatterplot(frbs, x='mag', y='P_Ox', s=3)

    # High confidednce
    mag_lim = 20.
    PATH_lim = 0.9
    high = frbs.P_Ox > PATH_lim
    bright = frbs.mag < mag_lim

    high_conf = frbs[high & bright]

    ax.text(0.05, 0.05, f'N(P_Ox>{PATH_lim}; m<{mag_lim})={len(high_conf)}/{len(frbs)}', 
            fontsize=14, transform=ax.transAxes)


    # Finish
    ffutils.set_fontsize(ax, 16.)

    pad = 0.2
    plt.tight_layout(pad=0.2,h_pad=pad,w_pad=pad)

    plt.savefig(outfile, dpi=300)#, bbox_inches='tight')
    plt.close()
    print(f"Wrote: {outfile}")

def fig_cosmos_ex(frb_file:str, outfile:str, local:str,
                  imsize:float=60./3600): 
    radec_sigma = [float(x)/3600. for x in local.split('x')]
    # Load COSMOS
    cosmos_df = cosmos.load_galaxies()

    # Load FRBs
    frbs = pandas.read_csv(frb_file)

    # Figures (4 subplots)
    fig = plt.figure(figsize=(8, 8))
    gs = gridspec.GridSpec(2,2)

    #for ss in range(4):
    for ss in range(3):
        ax = plt.subplot(gs[ss])

        # Select
        if ss == 0:
            # Bright one
            idx = np.argmin(np.abs(frbs.mag-19.))
        elif ss == 1:
            # Faint one
            idx = np.argmin(np.abs(frbs.mag-22.5))
        elif ss == 2:
            # Large offset
            idx = np.argmin(np.abs(frbs.loc_off-25.))

        # Continue
        frb = frbs.iloc[idx]
        center = SkyCoord(ra=frb.true_ra, dec=frb.true_dec, unit='deg')
        obs_center = SkyCoord(ra=frb.ra, dec=frb.dec, unit='deg')

        # Grab the galaxies
        gd_ra = cosmos_df['ra'].between(center.ra.deg-imsize/2., center.ra.deg+imsize/2.)
        gd_dec = cosmos_df['dec'].between(center.dec.deg-imsize/2., center.dec.deg+imsize/2.)
        gd_gal = gd_ra & gd_dec
        gal = cosmos_df[gd_gal]

        sns.scatterplot(gal, x='ra', y='dec', s=4, marker='o', ax=ax)

        ax.plot(center.ra.deg, center.dec.deg, 'rs', ms=3,
                fillstyle='none')

        # FRB and loclization
        ax.plot(obs_center.ra.deg, obs_center.dec.deg, 'k*', ms=2,
                fillstyle='none')
        # Ellipse
        ellipse = Ellipse(xy=(obs_center.ra.deg, obs_center.dec.deg),
                          width=2*radec_sigma[0], height=2*radec_sigma[1],
                          angle=0, color='k', fill=False)
        ax.add_artist(ellipse)

        # Fuss with the axes
        # Equalize
        ax.set_aspect('equal')

    pad = 0.2
    plt.tight_layout(pad=0.2,h_pad=pad,w_pad=pad)

    plt.savefig(outfile, dpi=300)#, bbox_inches='tight')
    plt.close()
    print(f"Wrote: {outfile}")

#### ########################## #########################
def main(pargs):

    if pargs.figure == 'scatter':
        # Scatter plot
        fig_scatter_Pm(f'PATH_{pargs.local}.csv', 
                     f'frb_monte_carlo_{pargs.local}.csv',
                     f'POx_mag_{pargs.local}.png')
    elif pargs.figure == 'roc':
        fig_roc(f'PATH_{pargs.local}.csv', 
                     f'frb_monte_carlo_{pargs.local}.csv',
                     f'ROC_{pargs.local}.png')
    elif pargs.figure == 'cosmos_ex':
        fig_cosmos_ex(f'frb_monte_carlo_{pargs.local}.csv',
                     f'fig_cosmos_ex_{pargs.local}.png',
                     pargs.local)


def parse_option():
    """
    This is a function used to parse the arguments in the training.
    
    Returns:
        args: (dict) dictionary of the arguments.
    """
    parser = argparse.ArgumentParser("CHIME/FRB GBO Figures")
    parser.add_argument("figure", type=str, 
                        help="function to execute: 'separation'")
    parser.add_argument('--local', type=str, default='3x15', help="Precitions (3x15)")
    parser.add_argument('--cmap', type=str, help="Color map")
    parser.add_argument('--debug', default=False, action='store_true',
                        help='Debug?')
    args = parser.parse_args()
    
    return args


# Command line execution
if __name__ == '__main__':

    pargs = parse_option()
    main(pargs)

# Examples
# python py/figs_chime_gbo.py cosmos_ex
# python py/figs_chime_gbo.py cosmos_ex --local 1x15

# Scatter
# python py/figs_chime_gbo.py scatter

# ROC
# python py/figs_chime_gbo.py roc