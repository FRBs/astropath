""" Figures related to CHIME/GBO calculations """
# imports
import numpy as np

import pandas
import argparse
from scipy import interpolate

import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Ellipse
from matplotlib.ticker import MultipleLocator

from astropy.coordinates import SkyCoord

from frb.figures import utils as ffutils
from astropath import cosmos

import analysis
import monte_carlo

from IPython import embed

def fig_mr_cdf(outfile:str='fig_mr_cdf.png', 
               chime_mr_tbl:str='CHIME_mr_5Jyms_150.parquet'): 
    """ Figure to check the m_r CDF

    Args:
        outfile (str, optional): Outfile. Defaults to 'fig_mr_cdf.png'.
        chime_mr_tbl (str, optional): Table of m_r CDF
        Defaults to 'CHIME_mr_5Jyms_150.parquet' generated by Mohit Bhardwaj (CHIME). 
    """
    # Load mr
    chime_mr = pandas.read_parquet(chime_mr_tbl)
    m_r_cuts = (chime_mr.m_r < monte_carlo.def_m_r_max) & (
            chime_mr.m_r > monte_carlo.def_m_r_min)
    chime_mr = chime_mr[m_r_cuts]

    # Generate a CDF
    n_chime = len(chime_mr)
    srt = np.argsort(chime_mr.m_r.values)

    f_cdf = interpolate.interp1d(
        np.arange(n_chime)/(n_chime-1),
        chime_mr.m_r.values[srt])

    # Figure
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111)

    # Plot cdf
    x_cdf = np.linspace(0., 1., 10000)
    y_cdf = f_cdf(x_cdf)
    ax.plot(y_cdf, x_cdf, 'k-', label='Full CDF')

    # Axes
    ax.set_xlabel(r'$m_r$')
    ax.set_ylabel('CDF')

    ffutils.set_fontsize(ax, 17.)
   
    pad = 0.2
    plt.tight_layout(pad=0.2,h_pad=pad,w_pad=pad)

    plt.savefig(outfile, dpi=300)#, bbox_inches='tight')
    plt.close()
    print(f"Wrote: {outfile}")


def fig_false_pos(path_file:str, frb_file:str,
                     outfile:str, local:str,
                     POX_cuts:float=[0.90, 0.95]): 
    """ False positive curve

    Args:
        path_file (str): PATH file
        frb_file (str): FRBs file
        outfile (str): Output file
        local (str): Description of localization
        POX_cut (float, optional): 
            P_Ox cut to use for the false positive curve
    """
    # parse
    frbs = analysis.parse_PATH(path_file, frb_file)

    def gen_cdf(tbl):
        
        # Sort on mag
        mag_srt = np.argsort(tbl.best_mag.values)

        # CDF
        cdf = np.cumsum(
            tbl.success.values[mag_srt]
            )/np.cumsum(np.ones(len(tbl)))

        cut_90 = np.argmin(np.abs(cdf - 0.9))
        cut_99 = np.argmin(np.abs(cdf - 0.99))
        print(f'90% completeness at {tbl.best_mag.values[mag_srt][cut_90]:g}')
        print(f'99% completeness at {tbl.best_mag.values[mag_srt][cut_99]:g}')

        return mag_srt, cdf 

    # Figure
    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111)

    # All
    mag_srt, cdf_all = gen_cdf(frbs)
    ax.plot(frbs.best_mag.values[mag_srt], 1.-cdf_all, 'g-', lw=2, label='All')

    # High P_Ox
    for POX_cut in POX_cuts:
        POx_cut = frbs.P_Ox > POX_cut
        mag_cut, cdf_cut = gen_cdf(frbs[POx_cut])
        ax.plot(frbs[POx_cut].best_mag.values[mag_cut], 1.-cdf_cut, 
                '-', lw=2, label=f'P(O|x) > {POX_cut*100}%')

    # Label
    ax.set_xlabel(r'$m_r$: Best Candidate')
    ax.set_ylabel(r'Cumulative FP Fraction ($m<m_r$)')
    ax.set_ylim(0., 1.)

    ax.legend()
    #ax.text(0.95, 0.05, f'Localization: {local}',
    #        fontsize=14, transform=ax.transAxes, ha='right')
    loc = MultipleLocator(base=0.1)
    #ax.xaxis.set_major_locator(loc)
    ax.yaxis.set_major_locator(loc)

    minorLocator = MultipleLocator(0.05)
    ax.yaxis.set_minor_locator(minorLocator)

    ax.grid(which='minor')
    ax.grid(which='major')

    ax.set_title(f'Localization: {local}', fontsize=16)

    # Finish
    ffutils.set_fontsize(ax, 16.)

    pad = 0.2
    plt.tight_layout(pad=0.2,h_pad=pad,w_pad=pad)

    plt.savefig(outfile, dpi=300)#, bbox_inches='tight')
    plt.close()
    print(f"Wrote: {outfile}")

    # Final stat
    #gd_mag = frbs.mag < 20.5
    gd_mag = frbs.mag < 22.0
    smpl = gd_mag & POx_cut
    print(f'We might target {np.sum(smpl)}/{len(frbs)}')

def fig_N_hosts(path_file:str, frb_file:str,
                     outfile:str, local:str,
                     POX_cuts:float=[0.85, 0.90, 0.925, 0.95, 0.97, 0.98, 0.99]): 
    """ Number of hosts (TP, FP) vs. P(O|x) cut

    Args:
        path_file (str): PATH file
        frb_file (str): FRBs file
        outfile (str): Output file
        local (str): Description of localization
        POX_cuts (float, optional): 
            P_Ox cuts 
    """
    # parse
    frbs = analysis.parse_PATH(path_file, frb_file)

    # Calculate
    TPs, FPs = [], []
    
    for POX_cut in POX_cuts:
        POx_cut = frbs.P_Ox > POX_cut
        # 
        TPs.append(np.sum(frbs[POx_cut].success.values))
        FPs.append(np.sum(np.invert(frbs[POx_cut].success.values)))

    # Table me
    df = pandas.DataFrame()
    df['P_Ox'] = POX_cuts
    df['N_TP'] = TPs
    df['N_FP'] = FPs


    # Figure
    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111)

    # Plot TPs
    ax.plot(df.P_Ox, df.N_TP/10000., 'o', color='blue', label='TP')
    ax.plot(df.P_Ox, df.N_FP/10000., 'x', color='red', label='FP')

    # Label
    ax.set_xlabel(r'$P(O|x)$ cut')
    ax.set_ylabel(r'Host Association: % of all KKO FRBs')
    ax.set_ylim(1./10000., max(df.N_TP)*1.3/10000.)
    ax.set_yscale('log')

    ax.legend(loc='lower left')
    '''
    loc = MultipleLocator(base=0.1)
    #ax.xaxis.set_major_locator(loc)
    ax.yaxis.set_major_locator(loc)

    minorLocator = MultipleLocator(0.05)
    ax.yaxis.set_minor_locator(minorLocator)

    ax.grid(which='minor')
    ax.grid(which='major')
    '''

    ax.set_title(f'Localization: {local}', fontsize=16)

    # Finish
    ffutils.set_fontsize(ax, 16.)

    pad = 0.2
    plt.tight_layout(pad=0.2,h_pad=pad,w_pad=pad)

    plt.savefig(outfile, dpi=300)#, bbox_inches='tight')
    plt.close()
    print(f"Wrote: {outfile}")


def fig_scatter_Pm(path_file:str, frb_file:str,
                     outfile:str): 
    """ Scatter plot of P_Ox vs. m_r

    Args:
        path_file (str): PATH file
        frb_file (str): FRBs file
        outfile (str): Output file
    """
    # parse
    frbs = analysis.parse_PATH(path_file, frb_file)

    #sns.set_style("whitegrid")

    fig = plt.figure(figsize=(12, 6))
    gs = gridspec.GridSpec(1,2)

    # Actual
    ax_actual = plt.subplot(gs[0])

    sns.scatterplot(frbs, x='mag', y='P_Ox', s=3, ax=ax_actual,
                    color='g')
    #ax_actual.set_yscale('log')
    ax_actual.set_ylim(0, 1.)

    ax_actual.set_xlabel(r'$m_r$:  True Galaxy')
    ax_actual.set_ylabel(r'$P(O|x)$:   Best Candidate')



    # Best
    correct = np.array(['yes']*len(frbs))
    correct[np.invert(frbs.success.values)] = 'no'
    frbs['correct'] = correct

    ax_best = plt.subplot(gs[1])
    sns.scatterplot(frbs, x='best_mag', y='P_Ox', s=3, 
                    ax=ax_best, hue='correct', 
                    palette='Set1')

    ax_best.set_xlabel(r'$m_r$:  Best Candidate')
    ax_best.set_ylabel(r'$P(O|x)$:   Best Candidate')
    ax_best.set_ylim(0, 1.)

    # High confidednce
    mag_lim = 20.
    PATH_lim = 0.9
    high = frbs.P_Ox > PATH_lim
    bright = frbs.mag < mag_lim

    high_conf = frbs[high & bright]

    #ax_best.text(0.05, 0.05, f'N(P_Ox>{PATH_lim}; m<{mag_lim})={len(high_conf)}/{len(frbs)}', 
    #        fontsize=14, transform=ax_best.transAxes)

    # Finish
    for ax in [ax_actual, ax_best]:
        ffutils.set_fontsize(ax, 16.)
        majorLocator = MultipleLocator(0.2)
        ax.yaxis.set_major_locator(majorLocator)
        minorLocator = MultipleLocator(0.05)
        ax.yaxis.set_minor_locator(minorLocator)

        ax.grid(which='minor')
        ax.grid(which='major')

    pad = 0.4
    plt.tight_layout(pad=0.2,h_pad=pad,w_pad=pad)

    plt.savefig(outfile, dpi=300)#, bbox_inches='tight')
    plt.close()
    print(f"Wrote: {outfile}")

def fig_cosmos_ex(frb_file:str, outfile:str, local:str,
                  imsize:float=60./3600): 
    """ Examlples of random FRBs in the COSMOS field

    Args:
        frb_file (str): FRBs file
        outfile (str): Output file
        local (str): Description of localization
        imsize (float, optional): Image size. Defaults to 60./3600.
    """
    radec_sigma = [float(x)/3600. for x in local.split('x')]
    # Load COSMOS
    cosmos_df = cosmos.load_galaxies()

    # PS
    in_ps = cosmos_df.mag_best.values < 23.2
    cosmos_ps = cosmos_df[in_ps].copy()
    cosmos_notps = cosmos_df[~in_ps].copy()
    

    # Load FRBs
    frbs = pandas.read_csv(frb_file)

    # Figures (4 subplots)
    fig = plt.figure(figsize=(8, 8))
    gs = gridspec.GridSpec(2,2)

    for ss in range(4):
        ax = plt.subplot(gs[ss])

        # Select
        if ss == 1:
            # Bright one
            lbl = 'Bright (r=19)'
            idx = np.argmin(np.abs(frbs.mag-19.))
        elif ss == 2:
            # Faint one
            lbl = 'Faint (r=22.5)'
            idx = np.argmin(np.abs(frbs.mag-22.5))
        elif ss == 3:
            # Large offset
            lbl = 'Large offset'
            idx = np.argmin(np.abs(frbs.loc_off-25.1))
        elif ss == 0:
            # Small offset
            lbl = 'Small offset'
            idx = np.argmin(frbs.loc_off)

        # Continue
        frb = frbs.iloc[idx]
        center = SkyCoord(ra=frb.true_ra, dec=frb.true_dec, unit='deg')
        obs_center = SkyCoord(ra=frb.ra, dec=frb.dec, unit='deg')

        # Grab the galaxies
        gd_ra = cosmos_df['ra'].between(center.ra.deg-imsize/2., center.ra.deg+imsize/2.)
        gd_dec = cosmos_df['dec'].between(center.dec.deg-imsize/2., center.dec.deg+imsize/2.)
        gd_gal = gd_ra & gd_dec

        ps_gal = cosmos_ps[gd_gal]
        notps_gal = cosmos_notps[gd_gal]

        # Galaxies
        # PS
        sns.scatterplot(ps_gal, x='ra', y='dec', s=8, marker='o', ax=ax)
        sns.scatterplot(notps_gal, x='ra', y='dec', s=3, marker='o', ax=ax,
                        color='lightgray')

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

        # Label
        ax.text(0.05, 0.05, f'{lbl}',
                fontsize=14, transform=ax.transAxes, ha='left')

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
    elif pargs.figure == 'false_pos':
        fig_false_pos(f'PATH_{pargs.local}.csv', 
                     f'frb_monte_carlo_{pargs.local}.csv',
                     f'FP_{pargs.local}.png',
                     pargs.local)
    elif pargs.figure == 'cosmos_ex':
        fig_cosmos_ex(f'frb_monte_carlo_{pargs.local}.csv',
                     f'fig_cosmos_ex_{pargs.local}.png',
                     pargs.local)
    elif pargs.figure == 'mr_cdf':
        fig_mr_cdf()
    elif pargs.figure == 'Nhosts':
        fig_N_hosts(f'PATH_{pargs.local}.csv', 
                     f'frb_monte_carlo_{pargs.local}.csv',
                     f'fig_Nhosts_{pargs.local}.png',
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
    parser.add_argument('--local', type=str, default='3x15', help="Predictions (3x15, 1x15)")
    parser.add_argument('--cmap', type=str, help="Color map")
    parser.add_argument('--debug', default=False, action='store_true',
                        help='Debug?')
    args = parser.parse_args()
    
    return args


# Command line execution
if __name__ == '__main__':

    pargs = parse_option()
    main(pargs)

# m_r CDF
# python py/figs_chime_gbo.py mr_cdf

# Examples
# python py/figs_chime_gbo.py cosmos_ex
# python py/figs_chime_gbo.py cosmos_ex --local 1x15

# Scatter P vs. m
# python py/figs_chime_gbo.py scatter

# False positives
# python py/figs_chime_gbo.py false_pos
# python py/figs_chime_gbo.py false_pos --local 1x15
# python py/figs_chime_gbo.py false_pos --local 005x15

# Number of hosts
# python py/figs_chime_gbo.py Nhosts 
# python py/figs_chime_gbo.py Nhosts --local 005x15