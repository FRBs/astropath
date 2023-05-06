""" Figures related to CHIME/GBO calculations """
# imports
import numpy as np

import pandas

import seaborn as sns
from matplotlib import pyplot as plt

from frb.figures import utils as ffutils

import analysis


def fig_simple_stats(path_file:str, frb_file:str,
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

# Command line execution
if __name__ == '__main__':

    # Generate FRBs
    fig_simple_stats('PATH_15x3.csv', 
                     'frb_monte_carlo_15x3.csv',
                     'simple_stats_15x3.png')