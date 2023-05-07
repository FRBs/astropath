""" Figures related to CHIME/GBO calculations """
# imports
import numpy as np

import pandas
import argparse

import seaborn as sns
from matplotlib import pyplot as plt

from frb.figures import utils as ffutils

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


#### ########################## #########################
def main(pargs):

    if pargs.figure == 'scatter':
        # Scatter plot
        fig_scatter_Pm('PATH_15x3.csv', 
                     'frb_monte_carlo_15x3.csv',
                     'POx_mag_15x3.png')
    # Gallery
    elif pargs.figure == 'roc':
        fig_roc('PATH_15x3.csv', 
                     'frb_monte_carlo_15x3.csv',
                     'ROC_15x3.png')

def parse_option():
    """
    This is a function used to parse the arguments in the training.
    
    Returns:
        args: (dict) dictionary of the arguments.
    """
    parser = argparse.ArgumentParser("CHIME/FRB GBO Figures")
    parser.add_argument("figure", type=str, 
                        help="function to execute: 'separation'")
    parser.add_argument('--frb', type=str, help="FRB name")
    parser.add_argument('--cmap', type=str, help="Color map")
    parser.add_argument('--debug', default=False, action='store_true',
                        help='Debug?')
    args = parser.parse_args()
    
    return args


# Command line execution
if __name__ == '__main__':

    pargs = parse_option()
    main(pargs)

# Scatter
# python py/figs_chime_gbo.py scatter

# ROC
# python py/figs_chime_gbo.py roc