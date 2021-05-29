""" Test methods in bayesian.py """

import os
from pkg_resources import resource_filename

import numpy as np
import pandas
import healpy as hp

from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units

from astropath import bayesian
from astropath import localization

import pytest

remote_data = pytest.mark.skipif(os.getenv('PATH_DATA') is None,
                                 reason='test requires dev suite')

def test_wcs():
    # Load up the localization
    lfile = os.path.join(resource_filename('astropath', 'tests'), 'files',
        'mask_frb201123_localization.fits.gz')
    hdul = fits.open(lfile)
    data = hdul[0].data
    header = hdul[0].header
    wcs = WCS(header)

    # Normalize
    data /= np.sum(data)
                         
    localiz = dict(type='wcs',
                   wcs_data=data, 
                   wcs_WCS=wcs)
    assert localization.vette_localization(localiz)
    
    # Approx center
    in_region = np.where(data > 0.)
    coord = wcs.pixel_to_world(in_region[1], in_region[0])
    cent_ra = np.mean(coord.ra.deg)
    cent_dec = np.mean(coord.dec.deg)

    # Calculate L_wx
    box_hwidth = 60.
    step_size = 1.
    ngrid = int(np.round(2*box_hwidth / step_size))
    x = np.linspace(-box_hwidth, box_hwidth, ngrid)
    xcoord, ycoord = np.meshgrid(x,x)
    ra = cent_ra + xcoord/3600. / np.cos(cent_dec*units.deg).value
    dec = cent_dec + ycoord/3600.

    L_wx = localization.calc_LWx(ra, dec, localiz)