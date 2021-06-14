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

from astropath import bayesian, priors
from astropath import localization

import pytest

remote_data = pytest.mark.skipif(os.getenv('PATH_DATA') is None,
                                 reason='test requires dev suite')

def test_raw_prior():
    # Inverse
    mag=np.array([21.,22.])
    ang_size = np.ones_like(mag)
    raw_PO = priors.raw_prior_Oi('inverse', ang_size, mag=np.array([21.,22.]))
    assert np.all(np.isclose(raw_PO, np.array([2736.80588898, 1074.04504448])))


def test_renorm_priors():
    # U=0
    raw_Oi = np.array([0.1, 0.2, 0.5])
    renorm_Oi = priors.renorm_priors(raw_Oi, 0.)

    assert np.all(renorm_Oi == np.array([0.125, 0.25 , 0.625]))

    # U != 0
    renorm_Oi = priors.renorm_priors(raw_Oi, 0.1)
    assert np.all(renorm_Oi == np.array([0.1125, 0.225, 0.5625]))
