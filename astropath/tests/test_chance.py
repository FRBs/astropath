""" Test methods in chance.py """
import numpy as np

from astropath import chance

import pytest


def test_driver_mag():
    # float
    rmag = 20.
    Pchance = chance.driver_sigma(rmag)
    assert np.isclose(Pchance, 0.00013535846962662733)
    assert isinstance(Pchance ,float)

    # array
    rmag = np.array([20., 21.])
    Pchance = chance.driver_sigma(rmag)
    assert np.isclose(Pchance, 0.00013535846962662733)
    assert isinstance(Pchance ,float)



