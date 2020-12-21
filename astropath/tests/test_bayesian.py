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

