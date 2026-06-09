""" Test the accuracy of the step size for the 
localization and galaxy convolution. 
"""

# imports
import os

import numpy as np
import matplotlib.pyplot as plt
#import pandas

from astropy.coordinates import SkyCoord
from astropy.coordinates import offset_by
from astropy import units
#from astropy.io import fits

from astropath import path
from astropath import localization
from astropath import bayesian

theta_prior = dict(max=6., PDF='exp', scale=1.)


def set_fontsize(ax, fsz):
    """
    Set the fontsize throughout an Axis

    Args:
        ax (Matplotlib Axis):
        fsz (float): Font size

    Returns:

    """
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                 ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(fsz)

def init_localization(a, b):

    Path = path.PATH()
    frb_coord = SkyCoord('21h44m25.255s -40d54m00.10s', frame='icrs')
    eellipse = dict(a=a, b=b, theta=0.)

    Path.init_localization('eellipse', 
                           center_coord=frb_coord, 
                           eellipse=eellipse)

    return Path


def big_localization():

    # FRB:  10" circle
    frb_local = 20. # arcsec
    path = init_localization(frb_local, frb_local)
    box_hwidth = 5 * frb_local # arcsec

    # Galaxy
    gal_coord = path.localiz['center_coord'].directional_offset_by(
        0.*units.deg, 1.0*units.arcsec)

    # Loop on galaxy size
    gal_sizes = np.array([0.1, 0.5, 1.0])

    #step_sizes = np.array([0.01, 0.05, 0.1, 0.2, 0.5, 1.])
    step_sizes = np.array([0.05, 0.1, 0.2, 0.5, 1., 2.]) # fast

    raw_pxO = np.zeros((len(gal_sizes), len(step_sizes)))
    cor_pxO = np.zeros((len(gal_sizes), len(step_sizes)))

    for ss, gal_size in enumerate(gal_sizes):

        # Loop on step size
        for ii, step_size in enumerate(step_sizes):
            print(f'Galaxy size: {gal_size}, Step size: {step_size}')

            # Calculate
            L_wx, p_wOi, grid_p, p_xOis = bayesian.px_Oi_fixedgrid(
                box_hwidth, path.localiz, np.array([gal_coord]),
                np.array([gal_size]), theta_prior, step_size=step_size, 
                return_debug=True)

            # Save
            raw_pxO[ss, ii] = p_xOis
            cor_pxO[ss, ii] = p_xOis / np.sum(p_wOi) / step_size**2
            #assert (np.sum(L_wx) * step_size**2) > 0.98, f'L_wx: {np.sum(L_wx) * step_size**2}'

    # Normalize by the first step size
    for ss in range(len(gal_sizes)):
        cor_pxO[ss, :] /= cor_pxO[ss, 0]
        raw_pxO[ss, :] /= raw_pxO[ss, 0]

    # Plot me
    fig = plt.figure(figsize=(10, 5))
    ax = plt.gca()

    clrs = ['red', 'green', 'blue']
    for ss, gal_size in enumerate(gal_sizes):
        ax.plot(step_sizes, cor_pxO[ss, :], color=clrs[ss],
                label=f'Galaxy size: {gal_size}')
        ax.plot(step_sizes, raw_pxO[ss, :], ls=':', color=clrs[ss])

    ax.legend(fontsize=15.)
    ax.set_xlabel('Step size (arcsec)')
    ax.set_ylabel('Normalized p_xO')

    ax.set_title(f'Big Localization, {frb_local}"')
    ax.set_xscale('log')
    #ax.set_yscale('log')
    #ax.grid(True)
    ax.grid(which='major', linewidth=0.8, alpha=0.7)
    ax.grid(which='minor', linewidth=0.5, alpha=0.3)
    set_fontsize(ax, 15.)

    # Vertical line at FRB local / 20
    ax.axvline(x=frb_local / 20., color='black', 
               linewidth=0.8)

    ax.set_ylim(0.9, 1.1)

    plt.tight_layout()
    plt.savefig('big_localization.png', dpi=300)
    print(f'Saved to big_localization.png')

# Command line
if __name__ == '__main__':
    big_localization()