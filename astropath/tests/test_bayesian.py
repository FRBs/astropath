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
    theta_prior = dict(max=6, PDF='uniform', scale=1.)
    assert priors.vet_theta_prior(theta_prior)
    pw_Oi_u = np.sum(bayesian.pw_Oi(theta, phi, theta_prior)) * grid_spacing_arcsec**2
    assert np.isclose(pw_Oi_u, 1., atol=1e-4)

    # Core
    theta_prior = dict(max=6, PDF='core', scale=1.)
    assert priors.vet_theta_prior(theta_prior)
    pw_Oi_c = np.sum(bayesian.pw_Oi(theta, phi, theta_prior)) * grid_spacing_arcsec**2
    assert np.isclose(pw_Oi_c, 1., atol=1e-4)

    # Exponential
    theta_prior = dict(max=6, PDF='exp', scale=1.)
    assert priors.vet_theta_prior(theta_prior)
    pw_Oi_e = np.sum(bayesian.pw_Oi(theta, phi, theta_prior)) * grid_spacing_arcsec**2
    assert np.isclose(pw_Oi_e, 1., atol=1e-4)

def test_error_ellipse():
    # This follows the FRB example Notebook

    # Set up localization
    frb_coord = SkyCoord('21h44m25.255s -40d54m00.10s', frame='icrs')
    eellipse = dict(a=0.1, b=0.1, theta=0.)
    localiz = dict(type='eellipse', center_coord=frb_coord, eellipse=eellipse)
    assert localization.vet_localization(localiz)

    # Candidates
    cand_file = os.path.join(resource_filename('astropath', 'data'), 'frb_example', 'frb180924_candidates.csv')
    candidates = pandas.read_csv(cand_file, index_col=0)
    c_coords = SkyCoord(ra=candidates.ra, dec=candidates.dec, unit='deg')

    # Priors
    theta_prior = dict(PDF='exp', max=6., scale=1.)
    cand_prior = dict(P_O_method='inverse', 
              P_U=0., 
              name='Adopted')
    assert priors.vet_theta_prior(theta_prior)
    assert priors.vet_cand_prior(cand_prior, candidates)

    # Raw priors
    raw_prior_Oi = priors.raw_prior_Oi(
        cand_prior['P_O_method'], candidates.half_light.values,
        mag=candidates.VLT_FORS2_g.values) 
    candidates['P_O_raw'] = raw_prior_Oi
    # Normalize
    candidates['P_O'] = priors.renorm_priors(
        candidates.P_O_raw.values, cand_prior['P_U'])

    # P(x|O)
    p_xOi = bayesian.px_Oi_fixedgrid(30.,  # box radius for grid, in arcsec
                       localiz,
                       c_coords,
                       candidates.half_light.values,
                       theta_prior,
                       step_size=0.02)
    candidates['p_xO'] = p_xOi
    
    # p(x)
    p_x = np.sum(candidates.P_O * candidates.p_xO)

    # Posteriors
    P_Oix = candidates.P_O * p_xOi / p_x
    candidates['P_Ox'] = P_Oix

    # Test
    assert np.isclose(candidates['P_Ox'].max(), 0.98951951218604)

def test_PU():
    # FRB
    frb_coord = SkyCoord('05h31m58.7013s +33d08m52.5536s', frame='icrs')
    eellipse = dict(a=0.1, b=0.1, theta=0.)
    localiz = dict(type='eellipse', center_coord=frb_coord, eellipse=eellipse)
    # Galaxies
    cand_file = os.path.join(resource_filename('astropath', 'data'), 'frb_example', 'frb121102_candidates.csv')
    candidates = pandas.read_csv(cand_file, index_col=0)
    c_coords = SkyCoord(ra=candidates.ra, dec=candidates.dec, unit='deg')

    candidates['sep'] = frb_coord.separation(c_coords).to('arcsec').value

    # Priors
    theta_prior = dict(PDF='uniform', max=6., scale=1.)
    cand_prior = dict(P_O_method='identical',
                  P_U=0.10,  # P(U)
                  name='Conservative')
    assert priors.vet_theta_prior(theta_prior)
    assert priors.vet_cand_prior(cand_prior, candidates)
    box_radius = 30.  # arcsec

    # Raw priors
    raw_prior_Oi = priors.raw_prior_Oi(cand_prior['P_O_method'],
                                         candidates.half_light.values,
                                         mag=candidates.GMOS_N_i.values)
    candidates['P_O_raw'] = raw_prior_Oi

    # Renormalize
    candidates['P_O'] = priors.renorm_priors(candidates.P_O_raw.values, 
                                               cand_prior['P_U'])
    assert np.isclose(np.sum(candidates.P_O), 0.9)

    # p(x|O)
    p_xOi = bayesian.px_Oi_fixedgrid(box_radius,  # box radius for grid, in arcsec
                           localiz,
                           c_coords,
                           candidates.half_light.values,
                           theta_prior)
    candidates['p_xO'] = p_xOi

    # p(x|U)
    p_xU = bayesian.px_U(box_radius)
    # p(x)
    p_x = np.sum(candidates.P_O * candidates.p_xO) + p_xU * cand_prior['P_U']
    assert np.isclose(p_x, 0.010006174202053927)

    # P(U|x)
    P_Ux = cand_prior['P_U'] * p_xU / p_x
    assert np.isclose(P_Ux, 0.0027760637799086035)

