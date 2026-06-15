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

def single_case(outfile:str, a, b, galaxy_size, 
                gal_offset:float=3., step_sizes=None,
                gal_pa:float=0.):

    box_hwidth = max(5 * a, gal_offset + 6 * galaxy_size) # arcsec

    # FRB
    path = init_localization(a, b)

    # Galaxy
    gal_coord = path.localiz['center_coord'].directional_offset_by(
        gal_pa*units.deg, gal_offset*units.arcsec)

    # Step sizes
    if step_sizes is None:
        step_sizes = np.array([0.02, 0.05, 0.1, 0.2, 0.5, 1., 2.]) # fast

    raw_pxO = np.zeros((1, len(step_sizes)))
    cor_pxO = np.zeros((1, len(step_sizes)))

    # Prescription
    best_step = max(b, galaxy_size)/20.

    for ss in [0]:

        # Loop on step size
        for ii, step_size in enumerate(step_sizes):
            print(f'Step size: {step_size}')

            # Calculate
            L_wx, p_wOi, grid_p, p_xOis = bayesian.px_Oi_fixedgrid(
                box_hwidth, path.localiz, np.array([gal_coord]),
                np.array([galaxy_size]), theta_prior, step_size=step_size, 
                return_debug=True)

            # Save
            raw_pxO[ss, ii] = p_xOis

            # Correction
            if b > galaxy_size:
                cor_pxO[ss, ii] = p_xOis / np.sum(p_wOi) / step_size**2
            else:
                cor_pxO[ss, ii] = p_xOis / np.sum(L_wx) / step_size**2
            #assert (np.sum(L_wx) * step_size**2) > 0.98, f'L_wx: {np.sum(L_wx) * step_size**2}'

    # Normalize by the first step size
    for ss in [0]:
        cor_pxO[ss, :] /= cor_pxO[ss, 0]
        raw_pxO[ss, :] /= raw_pxO[ss, 0]

    title = f'Localization, a={a}", b={b}; offset by {gal_offset}"; PA={gal_pa}"'
    plot_results(outfile, raw_pxO, cor_pxO, 
                 np.array([galaxy_size]), step_sizes, best_step,
                 var_name='Galaxy size', title=title)

def big_galaxy(galaxy_size:float=20.):

    box_hwidth = 6 * galaxy_size # arcsec

    # FRB  
    frb_coord = SkyCoord('21h44m25.255s -40d54m00.10s', frame='icrs')

    # Galaxy
    gal_coord = frb_coord.directional_offset_by(
        0.*units.deg, 3.0*units.arcsec)

    # Loop on FRB sizes
    frb_sizes = np.array([0.1, 0.5, 1.0, 5.])

    #step_sizes = np.array([0.01, 0.05, 0.1, 0.2, 0.5, 1.])
    step_sizes = np.array([0.05, 0.1, 0.2, 0.5, 1., 2.]) # fast

    raw_pxO = np.zeros((len(frb_sizes), len(step_sizes)))
    cor_pxO = np.zeros((len(frb_sizes), len(step_sizes)))

    for ss, frb_size in enumerate(frb_sizes):

        # Loop on step size
        for ii, step_size in enumerate(step_sizes):
            print(f'FRB size: {frb_size}, Step size: {step_size}')

            # FRB
            path = init_localization(frb_size, frb_size)

            # Calculate
            L_wx, p_wOi, grid_p, p_xOis = bayesian.px_Oi_fixedgrid(
                box_hwidth, path.localiz, np.array([gal_coord]),
                np.array([galaxy_size]), theta_prior, 
                step_size=step_size, 
                return_debug=True)

            # Save
            raw_pxO[ss, ii] = p_xOis
            cor_pxO[ss, ii] = p_xOis / np.sum(L_wx) / step_size**2
            # Another correction
            #cor_pxO[ss, ii] = cor_pxO[ss,ii] / np.sum(p_wOi) / step_size**2
            
            #assert (np.sum(L_wx) * step_size**2) > 0.98, f'L_wx: {np.sum(L_wx) * step_size**2}'

    # Normalize by the first step size
    for ss in range(len(frb_sizes)):
        cor_pxO[ss, :] /= cor_pxO[ss, 0]
        raw_pxO[ss, :] /= raw_pxO[ss, 0]

    # Plot
    outfile = 'big_galaxy.png'
    title = f'Big Galaxy, {galaxy_size}"'
    plot_results(outfile, raw_pxO, cor_pxO, 
                 frb_sizes, step_sizes, galaxy_size/20,
                 var_name='FRB size', title=title)

def big_localization(frb_local:float=20.):

    # FRB  
    path = init_localization(frb_local, frb_local)
    box_hwidth = 5 * frb_local # arcsec

    # Galaxy
    gal_coord = path.localiz['center_coord'].directional_offset_by(
        0.*units.deg, 1.0*units.arcsec)

    # Loop on galaxy size
    gal_sizes = np.array([0.1, 0.5, 1.0, 5.])

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

    # Plot
    outfile = 'big_localization.png'
    title = f'Big Localization, {frb_local}"'
    plot_results(outfile, raw_pxO, cor_pxO, 
                 gal_sizes, step_sizes, frb_local/20,
                 var_name='Galaxy size', title=title)


def gbo_localization(a=20., b=0.2):

    # FRB  
    path = init_localization(a, b)
    box_hwidth = 5 * a # arcsec

    # Galaxy
    gal_coord = path.localiz['center_coord'].directional_offset_by(
        0.*units.deg, 3.0*units.arcsec)

    # Loop on galaxy size
    gal_sizes = np.array([0.1, 0.5, 1.0, 5.])

    step_sizes = np.array([0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.])
    #step_sizes = np.array([0.05, 0.1, 0.2, 0.5, 1., 2.]) # fast

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

    # Plot
    outfile = 'gbo_localization.png'
    title = f'GBO Localization, a={a}", b={b}"'
    plot_results(outfile, raw_pxO, cor_pxO, 
                 gal_sizes, step_sizes, 0.1, 
                 var_name='Galaxy size', title=title)

def plot_results(outfile:str, raw_pxO, cor_pxO, 
                 var_sizes, step_sizes, 
                 pref_step_size:float,
                 var_name:str='Galaxy size',
                 title:str=None):

    # Plot me
    fig = plt.figure(figsize=(10, 5))
    ax = plt.gca()

    clrs = ['red', 'green', 'blue', 'purple']
    for ss, var_size in enumerate(var_sizes):
        ax.plot(step_sizes, cor_pxO[ss, :], color=clrs[ss],
                label=f'{var_name}: {var_size}')
        ax.plot(step_sizes, raw_pxO[ss, :], ls=':', color=clrs[ss])

    ax.legend(fontsize=15.)
    ax.set_xlabel('Step size (arcsec)')
    ax.set_ylabel('Normalized p_xO')

    if title is not None:
        ax.set_title(title, fontsize=15.)
    ax.set_xscale('log')
    #ax.set_yscale('log')
    #ax.grid(True)
    ax.grid(which='major', linewidth=0.8, alpha=0.7)
    ax.grid(which='minor', linewidth=0.5, alpha=0.3)
    set_fontsize(ax, 15.)

    # Vertical line 
    ax.axvline(x=pref_step_size, color='black', 
               linewidth=0.8)

    ax.set_ylim(0.9, 1.1)

    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    print(f'Saved to {outfile}')

# Command line
if __name__ == '__main__':

    #big_localization()
    #big_galaxy()
    #gbo_localization()

    # ############
    # Individual cases
    #single_case('single_case.png', 
    #            5., 5., 2., gal_offset=7.)
    single_case('single_case.png', 
                20., 0.2, 2., gal_offset=2., gal_pa=90.)