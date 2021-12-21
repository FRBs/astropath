""" Data model and code for Priors """
from typing import IO
import numpy as np

import pandas

from astropath import chance
from astropath import utils

theta_dmodel = {
    'PDF': dict(dtype=(str),
                options=['exp', 'core', 'uniform'],
                help='PDF shape for assumed transient offsets'),
    'max': dict(dtype=(float, np.floating, np.integer, int),
                help='Maximum offset (probability equals zero beyond this!)'),
}

cand_dmodel = {
    'P_O_method': dict(dtype=(str),
                options=['inverse', 'inverse_half', 'inverse_half2', 'identical'],
                help='Method for prior assignment of detected candidates.'),
    'P_U': dict(dtype=(float, np.floating),
                help='Prior for an unseen host.'),
    'name': dict(dtype=(str),
                help='Label for this prior.'),
}



def raw_prior_Oi(method, ang_size, mag=None, filter='r'):
    """
    Calculate raw prior for galaxy candidates
    If Sigma(m) is required, we adopt the Driver et al. 2016 
    evaluation 

    Args:
        method (str):
            inverse :: Assign inverse to Sigma_m
            inverse_ang :: Assign inverse to Sigma_m * half_light
            inverse_ang2 :: Assign inverse to Sigma_m * half_light**2
            identical :: All the same
        ang_size (float or np.ndarray):
            Angular size of the galaxy in arcsec
            Only required for several methods
        mag (float or np.ndarray, optional):
            Magnitudes of the sources;  assumes r-band and corrected
            for Galactic extinction

    Returns:
        float or np.ndarray:

    """
    # Convenience
    if method not in ['identical']:
        if filter != 'r':
            raise IOError("Not ready for this.  Best to go with what you have that is closest to r-band")
        Sigma_m = chance.driver_sigma(mag)

    # Do it
    if method == 'inverse':
        return 1. / Sigma_m
    elif method == 'inverse_ang':
        return 1. / Sigma_m / ang_size
    elif method == 'inverse_ang2':
        return 1. / Sigma_m / ang_size**2
    elif method == 'identical':
        return np.ones_like(ang_size)
    else:
        raise IOError("Bad method {} for prior_Oi".format(method))


def renorm_priors(raw_Oi, U):
    """
    Simple method to normalize the Priors

    Args:
        raw_Oi (np.ndarray):
            Unnormalized P_O values
        U (float):
            Prior for the galaxy being undetected

    Returns:
        np.ndarray: Normalized priors

    """
    raw_sum = np.sum(raw_Oi)
    return (1.-U) * raw_Oi/raw_sum

def vet_cand_prior(cand_prior:dict, candidates:pandas.DataFrame):
    """Vet the candidate prior dict

    Args:
        cand_prior (dict): 
        candidates (pandas.DataFrame): table of candidates

    Returns:
        bool: True if things check out ok
    """
    chk, disallowed_keys, badtype_keys = utils.vet_data_model(
        cand_prior, cand_dmodel)

    # Extras
    if cand_prior['P_O_method'] not in ['identical']:
        if 'mag' not in candidates.keys():
            print("You must provide candidate magnitudes to use this method!")

    # Return
    return chk

def vet_theta_prior(theta_prior:dict):
    """Vet the theta_prior dict

    Args:
        theta_prior (dict): 

    Returns:
        bool: True if things check out ok
    """
    chk, disallowed_keys, badtype_keys = utils.vet_data_model(
        theta_prior, theta_dmodel)
    return chk

def load_std_priors(nhalf:float=10., theta_max:float=6.):
    """
    Load the standard set of priors as defined
    in the PATH paper

    Args:
        nhalf (float, optional): Number of half-light radii for analysis. Defaults to 10..
        theta_max (float, optional): Maximum offset in half-light radii.
            Defaults to 6..

    Returns:
        dict: contains the 2 prior options:  "conservative" and "adopted"
    """
    # Theta priors
    theta_u = dict(method='uniform', max=theta_max)
    #theta_c = dict(method='core', max=theta_max)
    theta_e = dict(method='exp', max=theta_max)

    # Define them
    conservative = dict(theta=theta_u, O='identical', U=0., 
                        name='Conservative', nhalf=nhalf)
    adopted = dict(theta=theta_e, O='inverse', U=0., 
                   name='Adopted', nhalf=nhalf)

    priors = {}
    priors['conservative'] = conservative
    priors['adopted'] = adopted

    # Return
    return priors