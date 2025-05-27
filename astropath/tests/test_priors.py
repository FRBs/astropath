""" Test methods in bayesian.py """

import os
import numpy as np
import pandas
import healpy as hp

from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS

from unittest.mock import patch

from astropath import bayesian, priors
from astropath import localization

import pytest

remote_data = pytest.mark.skipif(os.getenv('PATH_DATA') is None,
                                 reason='test requires dev suite')


# ================================
# DRIVER TESTS (default filter='r')
# ================================

@patch('astropath.chance.driver_sigma')
def test_raw_prior_inverse(mock_driver_sigma):
    mag = np.array([21., 22.])
    ang_size = np.ones_like(mag)
    mock_driver_sigma.return_value = np.array([0.0003655, 0.0009314])
    result = priors.raw_prior_Oi('inverse', ang_size, mag=mag)
    expected = 1. / mock_driver_sigma.return_value
    np.testing.assert_allclose(result, expected, rtol=1e-5)

@patch('astropath.chance.driver_sigma')
def test_raw_prior_inverse_ang(mock_driver_sigma):
    mag = np.array([21., 22.])
    ang_size = np.array([2.0, 4.0])

    mock_driver_sigma.return_value = np.array([0.0003655, 0.0009314])
    result = priors.raw_prior_Oi('inverse_ang', ang_size, mag=mag)
    expected = 1. / mock_driver_sigma.return_value / ang_size
    np.testing.assert_allclose(result, expected, rtol=1e-5)

@patch('astropath.chance.driver_sigma')
def test_raw_prior_inverse_ang2(mock_driver_sigma):
    mag = np.array([21., 22.])
    ang_size = np.array([2.0, 4.0])

    mock_driver_sigma.return_value = np.array([0.0003655, 0.0009314])
    result = priors.raw_prior_Oi('inverse_ang2', ang_size, mag=mag)
    expected = 1. / mock_driver_sigma.return_value / ang_size**2
    np.testing.assert_allclose(result, expected, rtol=1e-5)


# ================================
# WINDHORST TESTS (filter='F200W')
# ================================

@patch('astropath.chance.windhorst_sigma')
def test_raw_prior_windhorst(mock_windhorst_sigma):
    mag = np.array([20., 22.])
    ang_size = np.ones_like(mag)

    mock_windhorst_sigma.return_value = np.array([0.001, 0.002])
    result = priors.raw_prior_Oi('inverse', ang_size, mag=mag, filter='F200W')
    expected = 1. / mock_windhorst_sigma.return_value
    np.testing.assert_allclose(result, expected, rtol=1e-5)

@patch('astropath.chance.windhorst_sigma')
def test_raw_prior_windhorst_ang(mock_windhorst_sigma):
    mag = np.array([20., 22.])
    ang_size = np.array([2.0, 4.0])

    mock_windhorst_sigma.return_value = np.array([0.001, 0.002])
    result = priors.raw_prior_Oi('inverse_ang', ang_size, mag=mag, filter='F200W')
    expected = 1. / mock_windhorst_sigma.return_value / ang_size
    np.testing.assert_allclose(result, expected, rtol=1e-5)

@patch('astropath.chance.windhorst_sigma')
def test_raw_prior_windhorst_ang2(mock_windhorst_sigma):
    mag = np.array([20., 22.])
    ang_size = np.array([2.0, 4.0])

    mock_windhorst_sigma.return_value = np.array([0.001, 0.002])
    result = priors.raw_prior_Oi('inverse_ang2', ang_size, mag=mag, filter='F200W')
    expected = 1. / mock_windhorst_sigma.return_value / ang_size**2
    np.testing.assert_allclose(result, expected, rtol=1e-5)


# ================================
# GENERIC TESTS
# ================================

def test_raw_prior_identical():
    ang_size = np.array([1.0, 2.0, 3.0])
    result = priors.raw_prior_Oi('identical', ang_size)
    np.testing.assert_array_equal(result, np.ones_like(ang_size))

def test_renorm_priors():
    raw_Oi = np.array([0.1, 0.2, 0.5])
    renorm_Oi = priors.renorm_priors(raw_Oi, 0.)
    assert np.allclose(renorm_Oi, np.array([0.125, 0.25, 0.625]))

    renorm_Oi = priors.renorm_priors(raw_Oi, 0.1)
    assert np.allclose(renorm_Oi, np.array([0.1125, 0.225, 0.5625]))

def test_norm_theta_priors():
    ngrid = 2000
    x = np.linspace(-30., 30., ngrid)
    grid_spacing_arcsec = x[1]-x[0]

    xcoord, ycoord = np.meshgrid(x, x)
    theta = np.sqrt(xcoord**2 + ycoord**2)

    for pdf in ['core', 'uniform', 'exp']:
        theta_prior = dict(max=6, PDF=pdf, scale=1.)
        p_wOi = bayesian.pw_Oi(theta, 2.0, theta_prior)
        assert np.isclose(np.sum(p_wOi) * grid_spacing_arcsec**2, 1., atol=1e-4), f'Failed normalization on {pdf}'

