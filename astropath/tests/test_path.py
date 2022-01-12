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

def test_frb():
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
    Path.init_theta_prior('exp', 6., 1.)

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


def test_gw():
    Path = path.PATH()
    # Load up the healpix
    lfile = os.path.join(resource_filename('astropath', 'data'), 'gw_examples',
                         'GW170817_skymap.fits.gz')
    gw170817 = hp.read_map(lfile)
    header = fits.open(lfile)[1].header

    Path.init_localization('healpix',
                   healpix_data=gw170817,
                   healpix_nside=header['NSIDE'],
                   healpix_ordering='NESTED',
                   healpix_coord='C')

    # Canidates
    galfile = os.path.join(resource_filename('astropath', 'data'), 'gw_examples',
                           'GW170817_galaxies.csv')
    cut_galaxies = pandas.read_csv(galfile)

    # Coords
    Path.init_candidates(cut_galaxies.RAJ2000, 
                         cut_galaxies.DEJ2000, 
                         cut_galaxies.maj.values * 60.,  # arcsec
                         mag=cut_galaxies.Bmag.values)

    # Candidate prior
    Path.init_cand_prior('inverse', P_U=0.)

    # Offset prior
    Path.init_theta_prior('exp', 6., 1.)

    # Priors
    p_O = Path.calc_priors()
    # Posterior
    P_Ox, P_Ux = Path.calc_posteriors('local', box_hwidth=30.)

    # Calculate p(x|O)
    assert np.isclose(np.max(Path.p_xOi), 0.05151238959823025)
    assert np.isclose(P_Ox.max(), 0.9999950450381409)