""" Defines localization datamodel for PATH anlysis and more"""
from typing import IO
import numpy as np

from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units

from IPython import embed

# Localization Data Model
localization_dmodel = {
    'type': dict(dtype=(str),
                options=['eellipse', 'WCS', 'healpix'],
                help='Localization type.'),
    'frb_coord': dict(dtype=(SkyCoord),
                help='FRB "central" coordinate'),
    'frb_eellipse': dict(dtype=(dict),
                help='FRB error ellipse with keys "a [arcsec]", "b [arcsec]", "theta [deg]"'),
    'healpix_data': dict(dtype=(np.ndarray, Table),
                help='Data containing the healpix information.'),
    'healpix_ordering': dict(dtype=(str),
                options=['NESTED', 'NUNIQ'],
                help='Ordering of the healpix information.'),
    'healpix_coord': dict(dtype=(str),
                options=['C'], # Celestial
                help='Coordinate system of the healpix.'),
    }

def calc_LWx(ra, dec, localization):
    if localization['type'] == 'eellipse':
        # Setup
        eellipse = localization['frb_eellipse']  # convenience
        pa_ee = eellipse['theta'] # PA of FRB error ellipse on the sky; deg
        dtheta = 90. - pa_ee  # Rotation to place the semi-major axis "a" of the ellipse along the x-axis we define
        #
        coord = SkyCoord(ra=ra, dec=dec, unit='deg')
        # Rotate to the FRB frame
        r = localization['frb_coord'].separation(coord).to('arcsec')
        pa_gal = localization['frb_coord'].position_angle(coord).to('deg')
        new_pa_gal = pa_gal + dtheta * units.deg
        # x, y of the box in FRB frame with x along major axis
        x_box = -r.value * np.sin(new_pa_gal).value
        y_box = r.value * np.cos(new_pa_gal).value

        # (orient semi-major axis "a" on our x axis)
        L_wx = np.exp(-x_box ** 2 / (2 * eellipse['a'] ** 2)) * np.exp(
            -y_box ** 2 / (2 * eellipse['b'] ** 2)) / (2*np.pi*eellipse['a']*eellipse['b'])
        #from pypeit.display import display
        #display.show_image(L_wx)
        #import pdb; pdb.set_trace()

    # Return
    return L_wx
        

def vette_localization(localization):
       
    chk = True
    # Loop on the keys
    disallowed_keys = []
    badtype_keys = []
    for key in localization.keys():
        # In data model?
        if not key in localization_dmodel.keys():
            disallowed_keys.append(key)
            chk = False
        # Check datat type
        if not isinstance(localization[key], 
                          localization_dmodel[key]['dtype']):
            badtype_keys.append(key)
            chk = False        

    # More
    if localization['type'] == 'eellipse':
        required_keys = ['frb_eellipse', 'frb_coord']
    elif localization['type'] == 'healpix':
        required_keys = ['healpix_data', 'healpix_ordering', 'healpix_coord']
    else:
        raise IOError("need required keys")
    for key in required_keys:
        if key not in localization.keys():
            chk = False
    
    # Return
    return chk