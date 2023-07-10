""" Definitions and methods related to COSMOS """

import os
from pkg_resources import resource_filename
import numpy as np

import pandas

from astropy import units
from astropy.coordinates import SkyCoord 

from IPython import embed

def cosmos_defs():
    """ Definitions for COSMOS """
    cdefs = dict(plate_scale=0.05,
                 filter='ACS_i',
                 )
    return cdefs

def load_galaxies(cosmos_file:str=None): 
    """ Load up the galaxies for COSMOS and set a few things

    The default input file may be found here:
    https://drive.google.com/file/d/1OoKa9ua_1XZJeH5NfCDBhEEFiwG2G2Bk/view?usp=sharing 
    

    Args:
        cosmos_file (str, optional): 
            COSMOS galaxy file. Defaults to 'cosmos_acs_iphot_200709.feather'.


    Returns:
        pandas.DataFrame: Table of galaxies
    """
    if cosmos_file is None:
        cosmos_file = os.path.join(resource_filename('astropath', 'data'), 'COSMOS', 
                                   'cosmos_acs_iphot_200709.feather')
        if not os.path.isfile(cosmos_file):
            raise IOError("You need to download the COSMOS galaxy file from GoogleDrive.  See the README.md file in that folder for info.")                                
    
    # Load
    df = pandas.read_feather(cosmos_file)

    # Cut
    df = df[df['mu_class'] == 1] # Filter out stars
    df = df[~np.isnan(df['kron_radius'])]

    pix_to_arcsec_kron = 0.03*units.arcsec
    df['half_light'] = df['kron_radius']*pix_to_arcsec_kron #arcsec

    # Return
    return df

