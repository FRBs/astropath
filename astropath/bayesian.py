"""Methods related to Bayesian association analysis"""

import numpy as np

import copy

from astropy import units
from astropy.coordinates import SkyCoord

from IPython import embed



def raw_prior_Oi(Pchance, Sigma_m, method):
    """
    Raw prior for a given set of Pchance values

    Proper normalization requires P(U) so that is done below

    Args:
        Pchance (np.ndarray):
            Chance probability
        Sigma_m (np.ndarray):
            Number density of sources on the sky brighter than m
        method (str):
            linear
            orig_inverse :: Assign inverse to P_chance
            inverse :: Assign inverse to Sigma_m
            identical :: All the same

    Returns:
        np.ndarray:

    """
    if method == 'linear':
        return 1 - Pchance
    elif method == 'orig_inverse':
        return 1. / Pchance
    elif method == 'inverse':
        return 1. / Sigma_m
    elif method == 'identical':
        return np.ones_like(Pchance)
    else:
        raise IOError("Bad method {} for prior_Oi".format(method))




def pw_Oi(r_w, theta_half, theta_prior, scale_half=1.):
    """
    Calculate p(w|O_i) for a given galaxy

    Args:
        r_w (np.ndarray):
            offset from galaxy center in arcsec
        theta_half (float):
            Half-light radius of this galaxy in arcsec
        theta_prior (dict):
            Parameters for theta prior
        scale_half (float):

    Returns:
        np.ndarray: Probability values; un-normalized

    """
    p = np.zeros_like(r_w)
    ok_w = r_w < theta_prior['max']*theta_half
    norm = 0
    if theta_prior['method'] == 'core':
        # Wolfram
        norm = theta_half * np.log(theta_prior['max']+1)
        if np.any(ok_w):
            p[ok_w] = theta_half / (r_w[ok_w] + theta_half) / norm
    elif theta_prior['method'] == 'uniform':
        norm = theta_half * theta_prior['max']
        if np.any(ok_w):
            p[ok_w] = 1. / norm
    elif theta_prior['method'] == 'exp':
        #norm = theta_half  - theta_half * (theta_prior['max']+1) * np.exp(-theta_prior['max'])
        norm = theta_half*scale_half * (scale_half - (
                scale_half+theta_prior['max'])*np.exp(-theta_prior['max']/scale_half))
        if np.any(ok_w):
            p[ok_w] = (r_w[ok_w] / theta_half) * np.exp(-r_w[ok_w]/(scale_half*theta_half)) / norm
    else:
        raise IOError("Bad theta method")
    #
    if norm == 0:
        raise ValueError("You forgot to normalize!")
    # Return
    return p


def px_Oi(box_radius, frb_coord, eellipse, cand_coords,
          theta_prior, step_size=0.1, return_grids=False):
    """
    Calculate p(x|O_i), the primary piece of the analysis

    Main concept:
        1. Set an area to analyze
        2. Discretize it to 0.1"
        3. Convolve

    Args:
        box_radius (float):
            Maximum radius for analysis, in arcsec
        eellipse (dict):
            Error ellipse for the FRB
            a, b, theta
        frb_coord (SkyCoord):
        cand_coords (SkyCoord):
            Coordinates of the candidate hosts
        theta_prior (dict):
            Parameters for theta prior
        step_size (float, optional):
            Step size for grid, in arcsec

    Returns:
        np.ndarray:

    """
    # Error ellipse
    pa_ee = eellipse['theta'] # deg
    dtheta = 90. - pa_ee  # Place a of ellipse along the x-axis
    # Set Equinox (for spherical offsets)
    frb_coord.equinox = cand_coords[0].equinox
    #
    ngrid = int(np.round(2*box_radius / step_size))
    x = np.linspace(-box_radius, box_radius, ngrid)
    xcoord, ycoord = np.meshgrid(x,x)

    # Build the grid around the FRB (orient semi-major axis "a" on our x axis)
    l_w = np.exp(-xcoord ** 2 / (2 * eellipse['a'] ** 2)) * np.exp(-ycoord ** 2 / (2 * eellipse['b'] ** 2))

    p_xMis, grids = [], []
    # TODO -- multiprocess this
    for icand, cand_coord in enumerate(cand_coords):

        # Calculate observed FRB location
        #dra, ddec = cand_coord.spherical_offsets_to(frb_coord)
        #xFRB = -dra.to('arcsec').value
        #yFRB = ddec.to('arcsec').value

        # #####################
        # l(w) -- 2D Gaussian

        # Rotate the galaxy
        r = frb_coord.separation(cand_coord).to('arcsec')
        pa_gal = frb_coord.position_angle(cand_coord).to('deg')
        new_pa_gal = pa_gal + dtheta * units.deg

        #r_wsq = (xcoord-xFRB)**2 + (ycoord-yFRB)**2
        #l_w = np.exp(-r_wsq/(2*sigR**2)) / sigR / np.sqrt(2*np.pi)

        # p(w|M_i)
        # x, y gal
        x_gal = -r.value * np.sin(new_pa_gal).value
        y_gal = r.value * np.cos(new_pa_gal).value
        r_w = np.sqrt((xcoord-x_gal)**2 + (ycoord-y_gal)**2)
        try:
            p_wMi = pw_Oi(r_w, theta_prior['r_half'][icand], theta_prior)
        except:
            embed(header='190 of bayesian')

        # Product
        grid_p = l_w * p_wMi

        # Save grids if returning
        if return_grids:
            grids.append(grid_p.copy())

        # Average
        p_xMis.append(np.mean(grid_p))

    # Return
    if return_grids:
        return np.array(p_xMis), grids
    else:
        return np.array(p_xMis)


def renorm_priors(raw_Oi, U):
    """
    Simple method to normalize the Priors

    Args:
        raw_Oi (np.ndarray):
        U (float):
            Prior for the FRB galaxy being undetected

    Returns:
        tuple: Normalized priors

    """
    raw_sum = np.sum(raw_Oi)
    return (1.-U) * raw_Oi/raw_sum


