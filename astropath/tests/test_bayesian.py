""" Test methods in bayesian.py """

import os
from pkg_resources import resource_filename

import numpy as np
import pandas
import healpy as hp

from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.io import fits

from astropath import bayesian
from astropath import localization

import pytest

remote_data = pytest.mark.skipif(os.getenv('PATH_DATA') is None,
                                 reason='test requires dev suite')

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

def test_error_ellipse():
    # Set up localization
    frb_coord = SkyCoord('21h44m25.255s -40d54m00.10s', frame='icrs')
    eellipse = dict(a=0.1, b=0.1, theta=0.)
    localiz = dict(type='eellipse', frb_coord=frb_coord, frb_eellipse=eellipse)
    assert localization.vette_localization(localiz)

    #


def test_PU():
    # FRB
    frb_coord = SkyCoord('05h31m58.7013s +33d08m52.5536s', frame='icrs')
    eellipse = dict(a=0.1, b=0.1, theta=0.)
    # Galaxies
    cand_file = os.path.join(resource_filename('astropath', 'data'), 'frb_example', 'frb121102_candidates.csv')
    candidates = pandas.read_csv(cand_file, index_col=0)
    c_coords = SkyCoord(ra=candidates.ra, dec=candidates.dec, unit='deg')

    candidates['sep'] = frb_coord.separation(c_coords).to('arcsec').value

    # Priors
    offset_prior = dict(method='uniform',
                        max=6.,
                        ang_size=candidates.half_light.values)
    priors = dict(offset=offset_prior,
                  O='identical',
                  U=0.10,  # P(U)
                  name='Conservative')
    box_radius = 30.  # arcsec

    # Raw priors
    raw_prior_Oi = bayesian.raw_prior_Oi(priors['O'],
                                         candidates.GMOS_N_i.values,
                                         half_light=candidates.half_light.values)
    candidates['P_O_raw'] = raw_prior_Oi

    # Renormalize
    candidates['P_O'] = bayesian.renorm_priors(candidates.P_O_raw.values, priors['U'])
    assert np.isclose(np.sum(candidates.P_O), 0.9)

    # p(x|O)
    p_xOi = bayesian.px_Oi(box_radius,  # box radius for grid, in arcsec
                           frb_coord,
                           eellipse,
                           c_coords,
                           priors['offset'])
    candidates['p_xO'] = p_xOi

    # p(x|U)
    p_xU = bayesian.px_U(box_radius)
    # p(x)
    p_x = np.sum(candidates.P_O * candidates.p_xO) + p_xU * priors['U']
    assert np.isclose(p_x, 0.010006174202053927)

    # P(U|x)
    P_Ux = priors['U'] * p_xU / p_x
    assert np.isclose(P_Ux, 0.0027760637799086035)


    pass


def test_healpix():
    # Load up the healpix
    lfile = os.path.join(resource_filename('astropath', 'data'), 'gw_examples',
                         'GW170817_skymap.fits.gz')
    gw170817 = hp.read_map(lfile)
    header = fits.open(lfile)[1].header

    # Galaxies
    galfile = os.path.join(resource_filename('astropath', 'data'), 'gw_examples',
                           'GW170817_galaxies.csv')
    cut_galaxies = pandas.read_csv(galfile)

    # Coords
    cut_gal_coord = SkyCoord(ra=cut_galaxies.RAJ2000, dec=cut_galaxies.DEJ2000, unit='deg')

    # PATH
    offset_prior = dict(method='exp',
                        max=6.,  # units of ang_size
                        ang_size=cut_galaxies.maj.values * 60.,  # arcsec
                        )
    priors = dict(offset=offset_prior,
                  O='inverse',
                  U=0.,
                  name='Adopted')

    # Priors
    raw_prior_Oi = bayesian.raw_prior_Oi(priors['O'],
                                         cut_galaxies.Bmag.values)
    cut_galaxies['P_O_raw'] = raw_prior_Oi
    cut_galaxies['P_O'] = bayesian.renorm_priors(cut_galaxies.P_O_raw.values, priors['U'])

    # Calculate p(x|O)
    p_xOi = bayesian.px_Oi_healpix(gw170817, header['NSIDE'],
                                   cut_gal_coord, offset_prior)  # , debug=True)
    np.isclose(np.max(p_xOi), 0.05150504596867476)



