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

from astropath.scripts import use_catalogs

import pytest

remote_data = pytest.mark.skipif(os.getenv('PATH_DATA') is None,
                                 reason='test requires dev suite')

def test_catalog():
    pargs = use_catalogs.parser(['128.680054,66.010750', '11.,11.,0.',
                                        '-U', '0.2', '--survey', 'Pan-STARRS'])
    Path = use_catalogs.main(pargs)

    # Test
    assert Path.candidates['P_Ux'][0] < 0.1