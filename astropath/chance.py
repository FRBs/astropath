""" Calculations related to P_chance """
import numpy as np

from scipy import interpolate

from importlib.resources import files
import os

from IPython import embed


# Spline parameters(globals) are for rmag vs sigma
driver_tck = (np.array([15., 15., 15., 15., 30., 30., 30., 30.]),
       np.array([-6.41580144, -3.53188049, -1.68500105, -0.63090954, 0., 0., 0., 0.]), 3)
driver_spl = interpolate.UnivariateSpline._from_tck(driver_tck)

# ---------------------------------------------------------------------------
# Updated Driver+16 counts: a spline fit directly to the Driver+16 Table 3 r-band
# DIFFERENTIAL counts (SDSS r, Subaru r, HST F606W)
# log10 of galaxies per mag per deg^2, valid 10 < r < 30. PATH's driver_spl above
# runs 1.3-2.1x high against Table 3 over 13 < r < 26.
# ---------------------------------------------------------------------------
driver16_tck = (np.array([10., 10., 10., 10., 30., 30., 30., 30.]),
              np.array([-2.24409004, 1.72855560, 5.13784270, 5.97024269, 0., 0., 0., 0.]), 3)
driver16_spl = interpolate.UnivariateSpline._from_tck(driver16_tck)
_D16_MMIN, _D16_MMAX = 10., 30.


def differential_driver16_sigma(rmag):
  """Galaxies per sq arcsec per mag at r = rmag, from the Driver+16 Table 3 fit.
  Clipped to the fit range 10 < r < 30."""
  m = np.clip(np.asarray(rmag, float), _D16_MMIN, _D16_MMAX)
  return 10**driver16_spl(m) / 3600.**2
  

def _build_driver16_cumulative(dm=1e-3):
  """Tabulate Sigma(<m) = integral of the differential counts, in sq arcsec^-1.
  Below r = 10 the counts continue as a power law matched to the log-slope of the
  fit at r = 10, so Sigma stays finite and positive for very bright galaxies."""
  m = np.arange(_D16_MMIN, _D16_MMAX + dm/2, dm)
  rho = differential_driver16_sigma(m) 
  k = (np.log10(rho[1]) - np.log10(rho[0])) / dm          # dex per mag at r = 10
  tail = rho[0] / (k * np.log(10.))                        # integral of the power law to r = 10
  cum = tail + np.concatenate([[0.], np.cumsum(0.5*(rho[1:] + rho[:-1]) * dm)])
  return m, cum, k, rho[0]
  
_D16_M, _D16_CUM, _D16_K, _D16_RHO10 = _build_driver16_cumulative()


def driver16_sigma(rmag):
  """Galaxies per sq arcsec with r < rmag (CUMULATIVE), from the Driver+16 Table 3 fit.
  This is the analogue of driver_sigma(), which the 'inverse' prior uses."""
  m = np.asarray(rmag, float)
  out = np.interp(m, _D16_M, _D16_CUM)
  bright = m < _D16_MMIN
  if np.any(bright):
      out = np.where(bright, _D16_RHO10 / (_D16_K*np.log(10.)) * 10**(_D16_K*(m - _D16_MMIN)), out)
  return out

def driver_sigma(mag):
    """
    Estimated incidence of galaxies per sq arcsec with r > rmag
    using Driver et al. 2016 number counts.

    Spline parameters (globals) are for rmag vs sigma

    Args:
        rmag (float or np.ndarray): r band magnitude of galaxy

    Returns:
        float or np.ndarray:  Galaxy number density

    """
    return 10**driver_spl(mag)

def windhorst_sigma(mag):
    """
    Estimated incidence of galaxies per sq arcsec with F200W > mag
    using Windhorst et al. 2024 number counts.

    Spline parameters (globals) are for F200W vs Num counts

    Args:
        mag (float or np.ndarray): F200W band magnitude of galaxy

    Returns:
        float or np.ndarray:  Galaxy number density

    """
    data_path = os.path.join(files('astropath'),'data','galaxy_num_counts','windhorst2023_F200W.npz')
    data = np.load(data_path)
    mag_f200w = data['mag']
    Num = data['Num(N/arcsec/0.5mag))']

    winhorst_spline = interpolate.interp1d(mag_f200w, Num, kind='cubic', fill_value='extrapolate')

    num_counts = winhorst_spline(mag)

    return num_counts


def differential_driver_sigma(rmag):
    """
    Estimated incidence of galaxies per sq arcsec per magnitude with r = rmag
    using Driver et al. 2016 number counts.

    Args:
        rmag (float or np.ndarray): r band magnitude of galaxy

    Returns:
        float or np.ndarray:  Galaxy differential number density rho

    """
    
    delta = 1e-5 # very small delta to numerically differentiate driver_spl
    rho = 10**driver_spl(rmag+delta) - 10**driver_spl(rmag-delta)
    
    # we calculate the density "per magnitude" by convention, hence must normalise
    # from 2 delta to unity
    rho /= 2.*delta
    
    return rho



def bloom_sigma(rmag):
    """
    Estimated incidence of galaxies per sq arcsec with r > rmag
    using Equation 3 of Bloom et al. 2002

    Args:
        rmag (float or np.ndarray): r band magnitude of galaxy

    Returns:
        float or np.ndarray:  Galaxy density

    """
    # Sigma(m)
    sigma = 1. / (3600. ** 2 * 0.334 * np.log(10)) * 10 ** (0.334 * (rmag - 22.963) + 4.320)
    return sigma


def pchance(rmag, sep, r_half, sigR, scale_rhalf=2., nsigma=2., ndens_eval='driver',
            orig_theff=False):
    """
    Calculate P_chance -- the probability of a chance association based
    on the number density of galaxies on the sky and the search area

    Args:
        rmag (np.ndarray):
            Extinction corrected apparent magnitude
        sep (np.ndarray):
            Angular separations; arcsec
        r_half (np.ndarray):
            Half light radii of the galaxies; arcsec
        sigR (float):
            1 sigma error in transient localization; assumed symmetric; arcsec
        scale_rhalf (float, optional):
            Weighting factor for the half-light radius
        nsigma:
            Weighting factor for the transient localization error
        ndens_eval (str, optinal):
            Number count source used
            'bloom': Hogg et al.
            'driver':  Driver et al. 2016
        orig_theff (bool, optional):
            Use the Bloom et al. 2002 prescription for theta_eff

    Returns:
        tuple: Pchance, nden -- chance probability and number density of sources given rmag

    """

    # Semi-original
    if orig_theff:
        Rs = np.stack([scale_rhalf * r_half, np.ones_like(r_half) * nsigma * sigR,
                       np.sqrt(sep ** 2 + (scale_rhalf * r_half) ** 2)])
        theta_eff = np.max(Rs, axis=0)
    else:
        # More conservative than usual
        theta_eff = np.sqrt(4*sigR**2 + 4*r_half**2 + sep**2)

    # Number density
    if ndens_eval =='bloom':
        nden = bloom_sigma(rmag)
    elif ndens_eval =='driver':
        nden = driver_sigma(rmag)
    else:
        raise IOError("Bad ndens evaluation")

    # Nbar
    Nbar = np.pi * theta_eff ** 2 * nden

    # Return Pchance and nden
    return 1. - np.exp(-Nbar), nden
