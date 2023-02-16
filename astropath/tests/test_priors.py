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

def test_norm_theta_priors():
    # Grid me
    ngrid = 2000
    x = np.linspace(-30., 30., ngrid)
    grid_spacing_arcsec = x[1]-x[0]

    xcoord, ycoord = np.meshgrid(x,x)
    theta = np.sqrt(xcoord**2 + ycoord**2)

    # Uniform
    for pdf in ['core', 'uniform', 'exp']:
        theta_prior = dict(max=6, PDF=pdf, scale=1.)
        p_wOi = bayesian.pw_Oi(theta,
                      2.0, # phi
                      theta_prior)
        assert np.isclose(np.sum(p_wOi)*grid_spacing_arcsec**2, 
                          1., atol=1e-4), f'Failed normalization on {pdf}'
    