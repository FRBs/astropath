"""Methods related to Bayesian association analysis"""

import numpy as np

from astropy import units
from astropy.coordinates import SkyCoord

from astropath import chance

from IPython import embed


def raw_prior_Oi(method, mag, half_light=None, Pchance=None):
    """
    Raw prior for a given set of magnitudes or Pchance values

    For the former, we adopt the Driver et al. 2016 evaluation of Sigma(m)

    Args:
        method (str):
            inverse :: Assign inverse to Sigma_m
            inverse1 :: Assign inverse to Sigma_m * half_light
            inverse2 :: Assign inverse to Sigma_m * half_light**2
            pchance_inverse :: Assign inverse to P_chance
            identical :: All the same
        mag (float or np.ndarray):
            Magnitudes of the sources;  assumed r-band
        half_light (float or np.ndarray, optional):
            Angular size of the galaxy in arcsec
            Only required for several methods
        Pchance (float or np.ndarray, optional):
            Chance probability
            Required for orig_inverse methods

    Returns:
        float or np.ndarray:

    """
    # Sigma_m -- not always used but that's ok
    Sigma_m = chance.driver_sigma(mag)

    # Do it
    if method == 'pchance_inverse':
        return 1. / Pchance
    elif method == 'inverse':
        return 1. / Sigma_m
    elif method == 'inverse1':
        if half_light is None:
            raise IOError("You must input angular size for method={}".format(method))
        return 1. / Sigma_m / half_light
    elif method == 'inverse2':
        if half_light is None:
            raise IOError("You must input angular size for method={}".format(method))
        return 1. / Sigma_m / half_light**2
    elif method == 'identical':
        return np.ones_like(Sigma_m)
    else:
        raise IOError("Bad method {} for prior_Oi".format(method))


def pw_Oi(r_w, phi, theta_prior, scale_half=1.):
    """
    Calculate p(w|O_i) for a given galaxy

    Args:
        r_w (np.ndarray):
            offset from galaxy center in arcsec
        phi (float):
            Angular size of the galaxy in arcsec
        theta_prior (dict):
            Parameters for theta prior
            Three methods are currently supported: core, uniform, exp
            See docs for further details
        scale_half (float, optional):
            Alternative scaling in the exponential

    Returns:
        np.ndarray: Probability values; un-normalized

    """
    p = np.zeros_like(r_w)
    ok_w = r_w < theta_prior['max']*phi
    if theta_prior['method'] == 'core':
        # Wolfram + Clancy
        norm = phi * np.log(theta_prior['max']/phi+1)
        if np.any(ok_w):
            p[ok_w] = phi / (r_w[ok_w] + phi) / norm
    elif theta_prior['method'] == 'uniform':
        norm = theta_prior['max']
        if np.any(ok_w):
            p[ok_w] = 1. / norm
    elif theta_prior['method'] == 'exp':
        norm = phi - np.exp(-scale_half*theta_prior['max']/phi) * (scale_half*theta_prior['max'] + phi)
        if np.any(ok_w):
            p[ok_w] = (r_w[ok_w] / phi) * np.exp(-scale_half*r_w[ok_w]/phi) / norm
    else:
        raise IOError("Bad theta method")
    #
    if norm == 0:
        raise ValueError("You forgot to normalize!")
    # Return
    return p


def px_Oi(box_hwidth, frb_coord, eellipse, cand_coords,
          theta_prior, step_size=0.1, return_grids=False):
    """
    Calculate p(x|O_i), the primary piece of the analysis

    Main concept:
        1. Set an area to analyze
        2. Discretize it to the step-size (e.g. 0.1")
        3. Convolve the FRB localization with the galaxy offset function

    Args:
        box_hwidth (float):
            Half-width of the analysis box, in arcsec
        frb_coord (SkyCoord):
        eellipse (dict):
            Error ellipse for the FRB
            a, b in arcsec, theta in deg
        cand_coords (SkyCoord):
            Coordinates of the candidate hosts
        theta_prior (dict):
            Parameters for theta prior
        step_size (float, optional):
            Step size for grid, in arcsec
        return_grids (bool, optional):
            if True, return the calcualtion grid

    Returns:
        np.ndarray or tuple: p(x|O_i) values and the grids if return_grids = True

    """
    # Error ellipse
    pa_ee = eellipse['theta'] # deg
    dtheta = 90. - pa_ee  # Place a of ellipse along the x-axis
    # Set Equinox (for spherical offsets)
    frb_coord.equinox = cand_coords[0].equinox
    #
    ngrid = int(np.round(2*box_hwidth / step_size))
    x = np.linspace(-box_hwidth, box_hwidth, ngrid)
    xcoord, ycoord = np.meshgrid(x,x)

    # Grid spacing
    grid_spacing_arcsec = x[1]-x[0]

    # #####################
    # Build the grid around the FRB (orient semi-major axis "a" on our x axis)
    # l(w) -- 2D Gaussian, normalized
    l_w = np.exp(-xcoord ** 2 / (2 * eellipse['a'] ** 2)) * np.exp(
        -ycoord ** 2 / (2 * eellipse['b'] ** 2)) / (2*np.pi*eellipse['a']*eellipse['b'])

    p_xOis, grids = [], []
    # TODO -- multiprocess this
    for icand, cand_coord in enumerate(cand_coords):


        # Rotate the galaxy
        r = frb_coord.separation(cand_coord).to('arcsec')
        pa_gal = frb_coord.position_angle(cand_coord).to('deg')
        new_pa_gal = pa_gal + dtheta * units.deg

        # p(w|O_i)
        # x, y gal
        x_gal = -r.value * np.sin(new_pa_gal).value
        y_gal = r.value * np.cos(new_pa_gal).value
        r_w = np.sqrt((xcoord-x_gal)**2 + (ycoord-y_gal)**2)
        p_wOi = pw_Oi(r_w, theta_prior['ang_size'][icand], theta_prior)

        # Product
        grid_p = l_w * p_wOi

        # Save grids if returning
        if return_grids:
            grids.append(grid_p.copy())

        # Sum
        p_xOis.append(np.sum(grid_p*grid_spacing_arcsec**2))

    # Return
    if return_grids:
        return np.array(p_xOis), grids
    else:
        return np.array(p_xOis)


def renorm_priors(raw_Oi, U):
    """
    Simple method to normalize the Priors

    Args:
        raw_Oi (np.ndarray):
        U (float):
            Prior for the FRB galaxy being undetected

    Returns:
        np.ndarray: Normalized priors

    """
    raw_sum = np.sum(raw_Oi)
    return (1.-U) * raw_Oi/raw_sum


