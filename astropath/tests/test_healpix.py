""" Test methods in bayesian.py """

import os
from pkg_resources import resource_filename

import numpy as np
import pandas
import healpy as hp

from astropy.table import Table
from astropy.coordinates import SkyCoord, Angle
from astropy.io import fits
from astropy import units

from astropath import healpix

import pytest

remote_data = pytest.mark.skipif(os.getenv('PATH_DATA') is None,
                                 reason='test requires dev suite')

def test_ellipse():
    frb_cen = SkyCoord(ra=326.1052292, dec=-40.9000278, unit='deg')
    pa_FRB = Angle(45, units.deg)
    a_FRB = Angle(0.4, units.arcsec)
    b_FRB = Angle(0.3, units.arcsec)

    # Run me
    hp_tbl = healpix.ellipse_to_healpix(frb_cen, pa_FRB, a_FRB, b_FRB)

    for key in ['UNIQ', 'PROBDENSITY']:
        assert key in hp_tbl.keys()

    assert len(hp_tbl) == 1783