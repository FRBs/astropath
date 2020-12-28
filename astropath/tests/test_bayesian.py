""" Test methods in bayesian.py """
import numpy as np

from astropath import bayesian

import pytest


def test_raw_prior():
    # Inverse
    raw_PO = bayesian.raw_prior_Oi('inverse', np.array([21.,22.]))
    assert np.all(np.isclose(raw_PO, np.array([2736.80588898, 1074.04504448])))


def test_renorm_priors():
    # U=0
    raw_Oi = np.array([0.1, 0.2, 0.5])
    renorm_Oi = bayesian.renorm_priors(raw_Oi, 0.)

    assert np.all(renorm_Oi == np.array([0.125, 0.25 , 0.625]))

    # U != 0
    renorm_Oi = bayesian.renorm_priors(raw_Oi, 0.1)
    assert np.all(renorm_Oi == np.array([0.1125, 0.225, 0.5625]))

def test_pw_Oi():
    # Generate a grid
    box_hwidth = 15. # arcsec
    step_size = 0.01 # arcsec
    phi = 1. # angular size in arcsec
    #
    ngrid = int(np.round(2*box_hwidth / step_size))
    x = np.linspace(-box_hwidth, box_hwidth, ngrid)
    xcoord, ycoord = np.meshgrid(x,x)

    grid_spacing_arcsec = x[1]-x[0]

    theta = np.sqrt(xcoord**2 + ycoord**2)

    # Uniform
    theta_prior = dict(max=6, method='uniform')
    pw_Oi_u = np.sum(bayesian.pw_Oi(theta, phi, theta_prior)) * grid_spacing_arcsec**2
    assert np.isclose(pw_Oi_u, 1., atol=1e-4)

    # Core
    theta_prior = dict(max=6, method='core')
    pw_Oi_c = np.sum(bayesian.pw_Oi(theta, phi, theta_prior)) * grid_spacing_arcsec**2
    assert np.isclose(pw_Oi_c, 1., atol=1e-4)

    # Exponential
    theta_prior = dict(max=6, method='exp')
    pw_Oi_e = np.sum(bayesian.pw_Oi(theta, phi, theta_prior)) * grid_spacing_arcsec**2
    assert np.isclose(pw_Oi_e, 1., atol=1e-4)

