""" Defines localization datamodel for PATH anlysis and more"""
from typing import IO
import numpy as np

import healpy as hp
import astropy_healpix

from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units
from astropy.wcs import WCS

from astropath import utils

from IPython import embed

# Localization Data Model
localization_dmodel = {
    'type': dict(dtype=(str),
                options=['eellipse', 'wcs', 'healpix'],
                help='Localization type.'),
    'center_coord': dict(dtype=(SkyCoord),
                help='"Central" coordinate'),
    'eellipse': dict(dtype=(dict),
                help='Error ellipse with keys "a [arcsec]", "b [arcsec]", "theta [deg]"'),
    'healpix_data': dict(dtype=(np.ndarray, Table),
                help='Data containing the healpix information.' \
            +' Input either as a simple numpy array for a full NESTED array' \
            +' or an astropy Table for NUNIQ format with' \
            +' columns UNIQ and PROBDENSITY.'),
    'healpix_nside': dict(dtype=(int),
                help='NSIDE value of healpix map.'),
    'healpix_ordering': dict(dtype=(str),
                options=['NESTED', 'NUNIQ'],
                help='Ordering of the healpix information.'),
    'healpix_coord': dict(dtype=(str),
                options=['C'], # Celestial
                help='Coordinate system of the healpix.'),
    'wcs_data': dict(dtype=(np.ndarray),
                help='PDF of localization.'),
    'wcs_WCS': dict(dtype=(WCS),
                help='WCS of the localization.'),
    }

def calc_LWx(ra:np.ndarray, dec:np.ndarray, localiz:dict):
    """Calculate the grid of L(w-x) values for
    the localization given an input ra, dec arrays

    Args:
        ra (np.ndarray): RA grid
        dec (np.ndarray): Dec grid
        localiz (dict): localization    
            Must follow the DataModel above

    Raises:
        IOError: [description]

    Returns:
        np.ndarray: L_wx grid
    """
            
    
    if localiz['type'] == 'eellipse':
        # Setup
        eellipse = localiz['eellipse']  # convenience
        pa_ee = eellipse['theta'] # PA of error ellipse on the sky; deg
        dtheta = 90. - pa_ee  # Rotation to place the semi-major axis "a" of the ellipse along the x-axis we define
        #
        coord = SkyCoord(ra=ra, dec=dec, unit='deg')
        coord.equinox = localiz['center_coord'].equinox
        # Rotate to the transient frame
        sep_box = localiz['center_coord'].separation(coord).to('arcsec')
        pa_box = localiz['center_coord'].position_angle(coord).to('deg')
        new_pa_box = pa_box + dtheta * units.deg
        # x, y of the box in transient frame with x along major axis
        x_box = -sep_box.value * np.sin(new_pa_box).value
        y_box = sep_box.value * np.cos(new_pa_box).value

        # Calculate
        L_wx = np.exp(-x_box ** 2 / (2 * eellipse['a'] ** 2)) * np.exp(
            -y_box ** 2 / (2 * eellipse['b'] ** 2)) / (2*np.pi*eellipse['a']*eellipse['b'])
    elif localiz['type'] == 'healpix':
        hp_index = hp.ang2pix(localiz['healpix_nside'], ra, dec, 
                              lonlat=True)
        # Healpix
        if localiz['healpix_ordering'] == 'NESTED':
            L_wx = localiz['healpix_data'][hp_index]
        elif localiz['healpix_ordering'] == 'NUNIQ':
            # Grab the pixels
            level, ipix = astropy_healpix.uniq_to_level_ipix(
                localiz['healpix_data']['UNIQ'])
            # Match
            match_id = utils.match_ids(hp_index.flatten(), ipix, require_in_match=False).reshape(ra.shape)
            gd_match = match_id > 0
            L_wx = np.zeros_like(ra)
            L_wx[gd_match] = localiz['healpix_data']['PROBDENSITY'][match_id[gd_match]]
            import pdb; pdb.set_trace()
        else:
            raise IOError("Bad healpix ordering! {}".format(localiz['healpix_ordering']))
    elif localiz['type'] == 'wcs':
        # Coords to col, row
        col_pix, row_pix = localiz['wcs_WCS'].all_world2pix(
            ra, dec, 1)
        # Recast
        col_pix = np.round(col_pix).astype(int)
        row_pix = np.round(row_pix).astype(int)
        # PDF
        L_wx = localiz['wcs_data'][row_pix, col_pix]
    else:
        raise IOError("Bad localization type!! {}".format(localiz['type']))

    # Return
    return L_wx
        

def vet_localization(localiz:dict):
    """Vet the localization object
    i.e. confirm it corresponds to the DataModel

    Args:
        localiz (dict): localization

    Raises:
        IOError: [description]

    Returns:
        bool: True if ok
    """

    chk, disallowed_keys, badtype_keys = utils.vet_data_model(localiz, localization_dmodel)

    # Required keys
    if localiz['type'] == 'eellipse':
        required_keys = ['eellipse', 'center_coord']
    elif localiz['type'] == 'healpix':
        required_keys = ['healpix_data', 'healpix_ordering', 
                         'healpix_coord', 'healpix_nside']
    elif localiz['type'] == 'wcs':
        required_keys = ['wcs_data', 'wcs_WCS']
    else:
        raise IOError("need required keys")
    for key in required_keys:
        if key not in localiz.keys():
            chk = False
    
    # Method specific tests
    if localiz['type'] == 'healpix':
        if localiz['healpix_ordering'] == 'NESTED':
            assert isinstance(localiz['healpix_data'], np.ndarray)
        elif localiz['healpix_ordering'] == 'NUNIQ':
            assert isinstance(localiz['healpix_data'], Table)
        
    # Return
    return chk