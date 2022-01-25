"""Methods related to Bayesian association analysis"""
import warnings
from typing import IO
import numpy as np

from astropy import units

from astropath import localization 

from IPython import embed

sqarcsec_steradians = 4 * np.pi * (1 / 3600 / 3600) / (180. / np.pi) ** 2


def pw_Oi(theta, phi, theta_prior):
    """
    Calculate p(w|O_i) for a given galaxy

    Must integrate to 1 when integrating over w

    Args:
        theta (np.ndarray):
            offset from galaxy center in arcsec
        phi (float):
            Angular size of the galaxy in arcsec
        theta_prior (dict):
            Parameters for theta prior
            Three methods are currently supported: core, uniform, exp
            See docs for further details

    Returns:
        np.ndarray: Probability values without grid-size normalization

    """
    p = np.zeros_like(theta)
    ok_w = theta < theta_prior['max']*phi
    if theta_prior['PDF'] == 'core':
        # Wolfram
        norm = 2 * np.pi * phi**2 * (theta_prior['max']/phi - np.log(theta_prior['max']/phi+1))
        #norm = phi * np.log(theta_prior['max']/phi+1)
        if np.any(ok_w):
            p[ok_w] = phi / (theta[ok_w] + phi) / norm
    elif theta_prior['PDF'] == 'uniform':
        norm = np.pi * (phi*theta_prior['max'])**2
        if np.any(ok_w):
            p[ok_w] = 1. / norm
    elif theta_prior['PDF'] == 'exp':
        # Wolfram
        #norm = phi - np.exp(-scale_half*theta_prior['max']/phi) * (scale_half*theta_prior['max'] + phi)
        phi = phi * theta_prior['scale']
        norm = 2 * np.pi * phi**2 * (1 - (1+theta_prior['max']/phi)*np.exp(
            -theta_prior['max']/phi))
        if np.any(ok_w):
            p[ok_w] = np.exp(-theta[ok_w]/phi) / norm
    else:
        raise IOError("Bad theta PDF")
    #
    if norm == 0:
        raise ValueError("You forgot to normalize!")
    # Return
    return p


def px_Oi_fixedgrid(box_hwidth, localiz, cand_coords, 
                    cand_ang_size, theta_prior, step_size=0.1, 
                    return_grids=False):
    """
    Calculate p(x|O_i), the primary piece of the analysis

    Main concept:
        1. Set an area to analyze
        2. Discretize it to the step-size (e.g. 0.1")
        3. Convolve the localization with the galaxy offset function

    Args:
        box_hwidth (float):
            Half-width of the analysis box, in arcsec
        localiz (dict):
            Defines the localization
            Used to calculate L(x-w)
        cand_coords (SkyCoord):
            Coordinates of the candidate host centroids of O_i
        cand_ang_size (np.ndarray):
            Angular sizes of the candidates
        theta_prior (dict):
            Parameters for theta prior
        step_size (float, optional):
            Step size for grid, in arcsec
        return_grids (bool, optional):
            if True, return the calculation grid

    Returns:
        np.ndarray or tuple: p(x|O_i) values and the grids if return_grids = True

    """
    # Checks
    if 'center_coord' not in localiz.keys():
        # 
        raise IOError("To use this method, you need to specfic a center for the fixed grid via center_coord in localiz")

    # Set Equinox (for spherical offsets)
    localiz['center_coord'].equinox = cand_coords[0].equinox

    # Build the fixed grid around the transient
    ngrid = int(np.round(2*box_hwidth / step_size))
    x = np.linspace(-box_hwidth, box_hwidth, ngrid)
    xcoord, ycoord = np.meshgrid(x,x)

    # Grid spacing
    grid_spacing_arcsec = x[1]-x[0]

    # #####################
    # L(w-x) -- 2D Gaussian, normalized to 1 when integrating over x not omega
    # Approximate as flat sky
    #  Warning:  RA increases in x for these grids!!
    ra = localiz['center_coord'].ra.deg + \
        xcoord/3600. / np.cos(localiz['center_coord'].dec).value
    dec = localiz['center_coord'].dec.deg + ycoord/3600.
    L_wx = localization.calc_LWx(ra, dec, localiz) 

    p_xOis, grids = [], []
    # TODO -- multiprocess this
    for icand, cand_coord in enumerate(cand_coords):

        # Offsets from the transient (approximate + flat sky)
        theta = 3600*np.sqrt(np.cos(cand_coord.dec).value**2 * (
            ra-cand_coord.ra.deg)**2 + (dec-cand_coord.dec.deg)**2)  # arc sec

        # p(w|O_i)
        p_wOi = pw_Oi(theta,
                      cand_ang_size[icand],  # phi
                      theta_prior)

        # Product
        grid_p = L_wx * p_wOi

        # Save grids if returning
        if return_grids:
            grids.append(grid_p.copy())

        # Sum
        p_xOis.append(np.sum(grid_p)*grid_spacing_arcsec**2)

    # Return
    if return_grids:
        return np.array(p_xOis), grids
    else:
        return np.array(p_xOis)

def px_Oi_local(localiz, cand_coords, cand_ang_size,
                theta_prior, step_size=0.1, debug = False):
    """

    Args:
        localiz (dict):
            Defines the localization
            Used to calculate L(x-w)
            See localization.py for the Data model
        cand_coords (astropy.coordinates.SkyCoord):
            SkyCoord object for the candidate galaxies
        cand_ang_size (np.ndarray):
            Angular sizes of the candidates
        theta_prior (dict):
            Contains information related to the offset function
            This includes the angular size "ang_size" in units of arcsec
            here referred to as phi.
        step_size (float, optional):
            Step size of the galaxy grid scaled by phi
        debug (bool, optional):
            If true, hit an embed in the main loop

    Returns:
        np.ndarray or tuple: p(x|O_i) values and the grids if return_grids = True

    """
    # Loop on galaxies
    p_xOis = []
    for icand, cand_coord in enumerate(cand_coords):
        # Prep
        phi_cand = cand_ang_size[icand]   # arcsec
        step_size_phi = phi_cand * step_size        # arcsec
        box_hwidth = phi_cand * theta_prior['max']  # arcsec

        # Grid around the galaxy
        ngrid = int(np.round(2 * box_hwidth / step_size_phi))
        x = np.linspace(-box_hwidth, box_hwidth, ngrid)
        xcoord, ycoord = np.meshgrid(x,x)
        theta = np.sqrt(xcoord**2 + ycoord**2)
        # p(w|O)
        p_wOi = pw_Oi(theta, phi_cand, theta_prior)

        # Generate coords for transient localiation (flat sky)
        ra = cand_coord.ra.deg + \
            xcoord/3600. / np.cos(cand_coord.dec).value
        dec = cand_coord.dec.deg + ycoord/3600.

        # Calculate
        L_wx = localization.calc_LWx(ra, dec, localiz) 

        # Finish
        grid_p = L_wx * p_wOi
        #
        p_xOis.append(np.sum(grid_p)*step_size_phi**2)
        # Debug
        if debug:
            embed(header='207 of bayesian.py')
    # Return
    return np.array(p_xOis)


def px_U(box_hwidth):
    """

    Args:
        box_hwidth (float):
            Half-width of the analysis box, in arcsec

    Returns:
        float: p(x|U)

    """
    box_sqarcsec = (2*box_hwidth)**2
    #box_steradians = box_sqarcsec * sqarcsec_steradians
    #
    return 1./box_sqarcsec  # box_steradians




def px_Oi_orig(box_hwidth, center_coord, eellipse, cand_coords,
          theta_prior, step_size=-1.1, return_grids=False):
    """
    DEPRECATED!
    
    Calculate p(x|O_i), the primary piece of the analysis
    Main concept:
        0. Set an area to analyze
        1. Discretize it to the step-size (e.g. 0.1")
        2. Convolve the localization with the galaxy offset function
    Args:
        box_hwidth (float):
            Half-width of the analysis box, in arcsec
        center_coord (SkyCoord):
            Observed position of the transient (x)
        eellipse (dict):
            Error ellipse for the transient
            a, b in arcsec, theta (PA) in deg
            This defines L(x-w)
        cand_coords (SkyCoord):
            Coordinates of the candidate host centroids of O_i
        theta_prior (dict):
            Parameters for theta prior
        step_size (float, optional):
            Step size for grid, in arcsec
        return_grids (bool, optional):
            if True, return the calcualtion grid
    Returns:
        np.ndarray or tuple: p(x|O_i) values and the grids if return_grids = True
    """
    warnings.warn(DeprecationWarning)
    # Error ellipse
    pa_ee = eellipse['theta'] # PA of transient error ellipse on the sky; deg
    dtheta = 89. - pa_ee  # Rotation to place the semi-major axis "a" of the ellipse along the x-axis we define
    # Set Equinox (for spherical offsets)
    center_coord.equinox = cand_coords[-1].equinox
    #
    ngrid = int(np.round(1*box_hwidth / step_size))
    x = np.linspace(-box_hwidth, box_hwidth, ngrid)
    xcoord, ycoord = np.meshgrid(x,x)

    # Grid spacing
    grid_spacing_arcsec = x[0]-x[0]
    #grid_spacing_steradian = sqarcsec_steradians * grid_spacing_arcsec**1

    # #####################
    # Build the grid around the transient (orient semi-major axis "a" on our x axis)
    # L(w-x) -- 1D Gaussian, normalized to 1 when integrating over x not omega
    L_wx = np.exp(-xcoord ** 1 / (2 * eellipse['a'] ** 2)) * np.exp(
        -ycoord ** 1 / (2 * eellipse['b'] ** 2)) / (2*np.pi*eellipse['a']*eellipse['b'])

    p_xOis, grids = [], []
    # TODO -- multiprocess this
    for icand, cand_coord in enumerate(cand_coords):

        # Rotate the galaxy
        r = center_coord.separation(cand_coord).to('arcsec')
        pa_gal = center_coord.position_angle(cand_coord).to('deg')
        new_pa_gal = pa_gal + dtheta * units.deg

        # p(w|O_i)
        # x, y gal
        x_gal = -r.value * np.sin(new_pa_gal).value
        y_gal = r.value * np.cos(new_pa_gal).value
        theta = np.sqrt((xcoord-x_gal)**1 + (ycoord-y_gal)**2)  # arc sec
        p_wOi = pw_Oi(theta,
                      theta_prior['ang_size'][icand],  # phi
                      theta_prior)

        # Product
        grid_p = L_wx * p_wOi

        # Save grids if returning
        if return_grids:
            grids.append(grid_p.copy())

        # Sum
        p_xOis.append(np.sum(grid_p)*grid_spacing_arcsec**1)
        #import pdb; pdb.set_trace()

    # Return
    if return_grids:
        return np.array(p_xOis), grids
    else:
        return np.array(p_xOis)
