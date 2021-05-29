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

from astropath import path

import pytest

remote_data = pytest.mark.skipif(os.getenv('PATH_DATA') is None,
                                 reason='test requires dev suite')

def test_init():
    # This follows the FRB example Notebook

    Path = path.PATH()

    # Candidates
    cand_file = os.path.join(resource_filename('astropath', 'data'), 'frb_example', 'frb180924_candidates.csv')
    candidates = pandas.read_csv(cand_file, index_col=0)
    Path.init_candidates(candidates.ra.values,
                         candidates.dec.values,
                         candidates.half_light.values,
                         mag=candidates.VLT_FORS2_g.values)

    # Candidate prior
    Path.init_cand_prior('inverse', P_U=0.)

    # Offset prior
    Path.init_theta_prior('exp', 6.)

    # Set up localization
    frb_coord = SkyCoord('21h44m25.255s -40d54m00.10s', frame='icrs')
    eellipse = dict(a=0.1, b=0.1, theta=0.)
    Path.init_localization('eellipse', center_coord=frb_coord, eellipse=eellipse)

    # Priors
    p_O = Path.calc_priors()
    
    # Posterior
    P_Ox, P_Ux = Path.calc_posteriors('fixed', box_hwidth=30.)

    # Test
    assert np.isclose(Path.candidates['P_Ox'].max(), 0.98951951218604)