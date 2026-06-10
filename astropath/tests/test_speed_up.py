""" Test the numpy speed-up of localization.calc_LWx (eellipse type).

The eellipse branch of localization.calc_LWx was re-implemented in pure
numpy (replacing astropy SkyCoord / separation / position_angle).  These
tests:

  1. Verify the new numpy code reproduces the original astropy code to
     within tight numerical tolerances.
  2. Time both implementations and report the speed-up, including a large
     7200x7200 pixel grid.

See prompts/speed_up.md for the analysis behind the change.
"""

import time

import numpy as np

from astropy.coordinates import SkyCoord
from astropy import units

from astropath import localization

import pytest


def _calc_LWx_eellipse_astropy(ra, dec, localiz):
    """Original (astropy-based) eellipse implementation of calc_LWx.

    Frozen here as the reference for regression/speed comparison.  This
    is a verbatim copy of the eellipse branch of localization.calc_LWx
    as it existed before the numpy rewrite.

    Args:
        ra (np.ndarray): RA grid (ICRS), deg.
        dec (np.ndarray): Dec grid (ICRS), deg.
        localiz (dict): Localization dict (eellipse type).

    Returns:
        np.ndarray: L(w-x) grid, same shape as ra/dec.
    """
    eellipse = localiz['eellipse']  # convenience
    pa_ee = eellipse['theta']  # PA of error ellipse on the sky; deg
    # Rotation to place the semi-major axis "a" along the x-axis
    dtheta = 90. - pa_ee
    #
    coord = SkyCoord(ra=ra, dec=dec, unit='deg')
    coord.equinox = localiz['center_coord'].equinox
    # Rotate to the transient frame
    sep_box = localiz['center_coord'].separation(coord).to('arcsec')
    pa_box = localiz['center_coord'].position_angle(coord).to('deg')
    new_pa_box = pa_box + dtheta * units.deg
    # x, y of the box in transient frame with x along major axis
    x_box = -sep_box.value * np.sin(new_pa_box).value
    y_box = sep_box.value * np.cos(new_pa_box).value
    # Calculate
    L_wx = np.exp(-x_box ** 2 / (2 * eellipse['a'] ** 2)) * np.exp(
        -y_box ** 2 / (2 * eellipse['b'] ** 2)) / (
        2 * np.pi * eellipse['a'] * eellipse['b'])
    return L_wx


def _build_grid(cent_ra, cent_dec, box_hwidth, ngrid):
    """Build an ra/dec grid around a center (flat-sky), as the callers do.

    Args:
        cent_ra (float): Central RA, deg.
        cent_dec (float): Central Dec, deg.
        box_hwidth (float): Half-width of the grid, arcsec.
        ngrid (int): Number of pixels per side.

    Returns:
        tuple: (ra, dec) numpy arrays of shape (ngrid, ngrid), deg.
    """
    x = np.linspace(-box_hwidth, box_hwidth, ngrid)
    xcoord, ycoord = np.meshgrid(x, x)
    # RA increases in x; flat-sky scaling by cos(dec)
    ra = cent_ra + xcoord / 3600. / np.cos(np.radians(cent_dec))
    dec = cent_dec + ycoord / 3600.
    return ra, dec


def _eellipse_localiz(cent_ra, cent_dec, a=1.0, b=0.5, theta=45.):
    """Construct an eellipse localization dict for testing.

    Args:
        cent_ra (float): Central RA, deg.
        cent_dec (float): Central Dec, deg.
        a (float): Semi-major axis, arcsec.
        b (float): Semi-minor axis, arcsec.
        theta (float): Position angle of the ellipse, deg (E of N).

    Returns:
        dict: Localization dict (eellipse type).
    """
    center = SkyCoord(ra=cent_ra, dec=cent_dec, unit='deg')
    return dict(type='eellipse',
                center_coord=center,
                eellipse=dict(a=a, b=b, theta=theta))


def test_eellipse_matches_astropy():
    """numpy calc_LWx must match the frozen astropy implementation."""
    cent_ra, cent_dec = 120.0, 32.0
    localiz = _eellipse_localiz(cent_ra, cent_dec)
    # Modest grid; box_hwidth=10", step_size=0.1" -> 200x200
    ra, dec = _build_grid(cent_ra, cent_dec, box_hwidth=10., ngrid=200)

    L_new = localization.calc_LWx(ra, dec, localiz)
    L_ref = _calc_LWx_eellipse_astropy(ra, dec, localiz)

    # Tight tolerance: the two should agree to roundoff
    assert np.allclose(L_new, L_ref, rtol=1e-8, atol=1e-12)
    # Sanity: L_wx peaks at the center pixel region and is finite
    assert np.all(np.isfinite(L_new))


def test_eellipse_matches_astropy_offcenter():
    """Match also at a southern, high-PA case to exercise the geometry."""
    cent_ra, cent_dec = 263.667, -50.768
    localiz = _eellipse_localiz(cent_ra, cent_dec, a=2.0, b=0.7, theta=110.)
    ra, dec = _build_grid(cent_ra, cent_dec, box_hwidth=15., ngrid=150)

    L_new = localization.calc_LWx(ra, dec, localiz)
    L_ref = _calc_LWx_eellipse_astropy(ra, dec, localiz)

    assert np.allclose(L_new, L_ref, rtol=1e-8, atol=1e-12)


def _time_call(func, ra, dec, localiz, reps=3):
    """Return the best-of-reps wall time (s) for func(ra, dec, localiz).

    Args:
        func (callable): Implementation to time.
        ra, dec (np.ndarray): Grid coords.
        localiz (dict): Localization dict.
        reps (int): Number of repetitions; the minimum is returned.

    Returns:
        float: Best wall-clock time in seconds.
    """
    best = np.inf
    for _ in range(reps):
        t0 = time.perf_counter()
        func(ra, dec, localiz)
        best = min(best, time.perf_counter() - t0)
    return best


@pytest.mark.parametrize("ngrid", [200, 1000, 7200])
def test_eellipse_speed_up(ngrid, capsys):
    """Time numpy vs astropy and report the speed-up.

    Includes the requested 7200x7200 pixel grid.  Also re-checks
    correctness on each grid (subsampled for the largest grid to keep
    memory/time reasonable).
    """
    cent_ra, cent_dec = 120.0, 32.0
    localiz = _eellipse_localiz(cent_ra, cent_dec)
    # Keep angular scale fixed-ish; box grows with ngrid
    box_hwidth = 0.05 * ngrid  # arcsec (0.1"/pix)
    ra, dec = _build_grid(cent_ra, cent_dec, box_hwidth, ngrid)

    t_np = _time_call(localization.calc_LWx, ra, dec, localiz)
    t_ap = _time_call(_calc_LWx_eellipse_astropy, ra, dec, localiz)

    # Correctness check (single eval)
    L_new = localization.calc_LWx(ra, dec, localiz)
    L_ref = _calc_LWx_eellipse_astropy(ra, dec, localiz)
    assert np.allclose(L_new, L_ref, rtol=1e-8, atol=1e-12)

    speedup = t_ap / t_np if t_np > 0 else np.inf
    # Report (use -s to see this output)
    with capsys.disabled():
        print(
            "\n  grid %5dx%-5d (%9d pts): "
            "astropy %8.2f ms  numpy %8.2f ms  speed-up %5.1fx"
            % (ngrid, ngrid, ngrid * ngrid,
               t_ap * 1e3, t_np * 1e3, speedup))

    # The numpy path should never be slower than astropy
    assert t_np <= t_ap
