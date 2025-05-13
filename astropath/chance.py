""" Calculations related to P_chance """
import numpy as np

from scipy import interpolate

from importlib.resources import files
import os

from IPython import embed


# Spline parameters(globals) are for rmag vs sigma
driver_tck = (np.array([15., 15., 15., 15., 30., 30., 30., 30.]),
       np.array([-6.41580144, -3.53188049, -1.68500105, -0.63090954, 0., 0., 0., 0.]), 3)
driver_spl = interpolate.UnivariateSpline._from_tck(driver_tck)



def driver_sigma(mag):
    """
    Estimated incidence of galaxies per sq arcsec with r > rmag
    using Driver et al. 2016 number counts.

    Spline parameters (globals) are for rmag vs sigma

    Args:
        rmag (float or np.ndarray): r band magnitude of galaxy

    Returns:
        float or np.ndarray:  Galaxy number density

    """
    return 10**driver_spl(mag)

def windhorst_sigma(mag):
    """
    Estimated incidence of galaxies per sq arcsec with F200W > mag
    using Windhorst et al. 2024 number counts.

    Spline parameters (globals) are for F200W vs Num counts

    Args:
        mag (float or np.ndarray): F200W band magnitude of galaxy

    Returns:
        float or np.ndarray:  Galaxy number density

    """
    data_path = os.path.join(files('astropath'),'data','galaxy_num_counts','windhorst2023_F200W.npz')
    data = np.load(data_path)
    mag_f200w = data['mag']
    Num = data['Num(N/arcsec/0.5mag))']

    winhorst_spline = interpolate.interp1d(mag_f200w, Num, kind='cubic', fill_value='extrapolate')

    num_counts = winhorst_spline(mag)

    return num_counts


def bloom_sigma(rmag):
    """
    Estimated incidence of galaxies per sq arcsec with r > rmag
    using Equation 3 of Bloom et al. 2002

    Args:
        rmag (float or np.ndarray): r band magnitude of galaxy

    Returns:
        float or np.ndarray:  Galaxy density

    """
    # Sigma(m)
    sigma = 1. / (3600. ** 2 * 0.334 * np.log(10)) * 10 ** (0.334 * (rmag - 22.963) + 4.320)
    return sigma


def pchance(rmag, sep, r_half, sigR, scale_rhalf=2., nsigma=2., ndens_eval='driver',
            orig_theff=False):
    """
    Calculate P_chance -- the probability of a chance association based
    on the number density of galaxies on the sky and the search area

    Args:
        rmag (np.ndarray):
            Extinction corrected apparent magnitude
        sep (np.ndarray):
            Angular separations; arcsec
        r_half (np.ndarray):
            Half light radii of the galaxies; arcsec
        sigR (float):
            1 sigma error in transient localization; assumed symmetric; arcsec
        scale_rhalf (float, optional):
            Weighting factor for the half-light radius
        nsigma:
            Weighting factor for the transient localization error
        ndens_eval (str, optinal):
            Number count source used
            'bloom': Hogg et al.
            'driver':  Driver et al. 2016
        orig_theff (bool, optional):
            Use the Bloom et al. 2002 prescription for theta_eff

    Returns:
        tuple: Pchance, nden -- chance probability and number density of sources given rmag

    """

    # Semi-original
    if orig_theff:
        Rs = np.stack([scale_rhalf * r_half, np.ones_like(r_half) * nsigma * sigR,
                       np.sqrt(sep ** 2 + (scale_rhalf * r_half) ** 2)])
        theta_eff = np.max(Rs, axis=0)
    else:
        # More conservative than usual
        theta_eff = np.sqrt(4*sigR**2 + 4*r_half**2 + sep**2)

    # Number density
    if ndens_eval =='bloom':
        nden = bloom_sigma(rmag)
    elif ndens_eval =='driver':
        nden = driver_sigma(rmag)
    else:
        raise IOError("Bad ndens evaluation")

    # Nbar
    Nbar = np.pi * theta_eff ** 2 * nden

    # Return Pchance and nden
    return 1. - np.exp(-Nbar), nden
