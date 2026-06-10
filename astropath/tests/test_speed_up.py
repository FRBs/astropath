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
from astropath import bayesian

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


# ---------------------------------------------------------------------------
# px_Oi_fixedgrid: replaced per-candidate / center astropy access with numpy
# (pre-extracted arrays).  Reference below freezes the original astropy
# implementation for regression + speed comparison.
# ---------------------------------------------------------------------------


def _px_Oi_fixedgrid_astropy(box_hwidth, localiz, cand_coords,
                             cand_ang_size, theta_prior, step_size=0.1):
    """Original (astropy-based) px_Oi_fixedgrid implementation.

    Verbatim copy of the loop/grid logic as it existed before the numpy
    rewrite (iterates the SkyCoord array, reads .ra/.dec per candidate,
    sets equinox).  Kept here as the reference for the tests.  Returns
    only the p(x|O_i) array (the only output the tests need).

    Args:
        box_hwidth (float): Half-width of the analysis box, arcsec.
        localiz (dict): Localization dict (must have center_coord).
        cand_coords (SkyCoord): Candidate host coordinates.
        cand_ang_size (np.ndarray): Candidate angular sizes, arcsec.
        theta_prior (dict): Offset-prior parameters.
        step_size (float): Grid step size, arcsec.

    Returns:
        np.ndarray: p(x|O_i) for each candidate.
    """
    # Set Equinox (for spherical offsets)
    localiz['center_coord'].equinox = cand_coords[0].equinox
    # Build the fixed grid around the transient
    ngrid = int(np.round(2 * box_hwidth / step_size))
    x = np.linspace(-box_hwidth, box_hwidth, ngrid)
    xcoord, ycoord = np.meshgrid(x, x)
    grid_spacing_arcsec = x[1] - x[0]
    # L(w-x); RA increases in x (flat-sky)
    ra = localiz['center_coord'].ra.deg + \
        xcoord / 3600. / np.cos(localiz['center_coord'].dec).value
    dec = localiz['center_coord'].dec.deg + ycoord / 3600.
    L_wx = localization.calc_LWx(ra, dec, localiz)
    p_xOis = []
    for icand, cand_coord in enumerate(cand_coords):
        # Offsets from the transient (flat sky)
        theta = 3600 * np.sqrt(np.cos(cand_coord.dec).value**2 * (
            ra - cand_coord.ra.deg)**2
            + (dec - cand_coord.dec.deg)**2)  # arcsec
        p_wOi = bayesian.pw_Oi(theta, cand_ang_size[icand], theta_prior)
        grid_p = L_wx * p_wOi
        p_xOis.append(np.sum(grid_p) * grid_spacing_arcsec**2)
    return np.array(p_xOis)


def _make_candidates(cent_ra, cent_dec, ncand=10, spread=8.):
    """Build a SkyCoord of candidates scattered around a center.

    Args:
        cent_ra (float): Central RA, deg.
        cent_dec (float): Central Dec, deg.
        ncand (int): Number of candidates.
        spread (float): Half-spread of the scatter, arcsec.

    Returns:
        tuple: (cand_coords SkyCoord, cand_ang_size ndarray arcsec).
    """
    # Deterministic offsets (no RNG) spanning +/- spread arcsec
    off = np.linspace(-spread, spread, ncand)
    cand_ra = cent_ra + off / 3600. / np.cos(np.radians(cent_dec))
    cand_dec = cent_dec + off[::-1] / 3600.
    cand_coords = SkyCoord(ra=cand_ra, dec=cand_dec, unit='deg')
    # Angular sizes 0.5-2.0 arcsec
    cand_ang_size = np.linspace(0.5, 2.0, ncand)
    return cand_coords, cand_ang_size


def test_px_Oi_fixedgrid_matches_astropy():
    """numpy px_Oi_fixedgrid must match the frozen astropy version."""
    cent_ra, cent_dec = 120.0, 32.0
    localiz = _eellipse_localiz(cent_ra, cent_dec, a=1.0, b=0.6, theta=30.)
    cand_coords, cand_ang_size = _make_candidates(cent_ra, cent_dec)
    theta_prior = dict(PDF='exp', max=6., scale=0.5)

    p_new = bayesian.px_Oi_fixedgrid(
        10., localiz, cand_coords, cand_ang_size, theta_prior)
    p_ref = _px_Oi_fixedgrid_astropy(
        10., localiz, cand_coords, cand_ang_size, theta_prior)

    # Same math, only coordinate extraction changed -> agree to roundoff
    assert np.allclose(p_new, p_ref, rtol=1e-10, atol=1e-15)


def test_px_Oi_fixedgrid_speed_up(capsys):
    """Time numpy vs astropy px_Oi_fixedgrid and report the speed-up."""
    cent_ra, cent_dec = 120.0, 32.0
    localiz = _eellipse_localiz(cent_ra, cent_dec, a=1.0, b=0.6, theta=30.)
    # Many candidates: the per-candidate astropy overhead is the target
    cand_coords, cand_ang_size = _make_candidates(
        cent_ra, cent_dec, ncand=50)
    theta_prior = dict(PDF='exp', max=6., scale=0.5)

    def _np(*_):
        return bayesian.px_Oi_fixedgrid(
            10., localiz, cand_coords, cand_ang_size, theta_prior)

    def _ap(*_):
        return _px_Oi_fixedgrid_astropy(
            10., localiz, cand_coords, cand_ang_size, theta_prior)

    # _time_call signature is (func, ra, dec, localiz); pass dummies
    t_np = _time_call(_np, None, None, None)
    t_ap = _time_call(_ap, None, None, None)

    assert np.allclose(_np(), _ap(), rtol=1e-10, atol=1e-15)

    speedup = t_ap / t_np if t_np > 0 else np.inf
    with capsys.disabled():
        print(
            "\n  px_Oi_fixedgrid (50 cand, 200x200): "
            "astropy %8.2f ms  numpy %8.2f ms  speed-up %5.1fx"
            % (t_ap * 1e3, t_np * 1e3, speedup))

    # numpy path should not be slower
    assert t_np <= t_ap


# ---------------------------------------------------------------------------
# Optional numba kernel (use_numba=True) for px_Oi_fixedgrid.
# ---------------------------------------------------------------------------

# Skip the numba tests entirely if numba is not installed
numba_required = pytest.mark.skipif(
    not bayesian.HAS_NUMBA, reason='numba not installed')


@numba_required
@pytest.mark.parametrize("pdf", ['exp', 'core', 'uniform'])
def test_px_Oi_fixedgrid_numba_matches_numpy(pdf):
    """use_numba=True must match the numpy path for every PDF."""
    cent_ra, cent_dec = 120.0, 32.0
    localiz = _eellipse_localiz(cent_ra, cent_dec, a=1.0, b=0.6, theta=30.)
    cand_coords, cand_ang_size = _make_candidates(cent_ra, cent_dec)
    theta_prior = dict(PDF=pdf, max=6., scale=0.5)

    p_np = bayesian.px_Oi_fixedgrid(
        10., localiz, cand_coords, cand_ang_size, theta_prior,
        use_numba=False)
    p_nb = bayesian.px_Oi_fixedgrid(
        10., localiz, cand_coords, cand_ang_size, theta_prior,
        use_numba=True)

    # Same math, just fused in the kernel -> agree to roundoff
    assert np.allclose(p_np, p_nb, rtol=1e-10, atol=1e-15)


def test_px_Oi_fixedgrid_numba_fallback_without_numba(monkeypatch):
    """use_numba=True falls back to numpy (with a warning) if numba is
    unavailable, and still returns the correct result."""
    cent_ra, cent_dec = 120.0, 32.0
    localiz = _eellipse_localiz(cent_ra, cent_dec, a=1.0, b=0.6, theta=30.)
    cand_coords, cand_ang_size = _make_candidates(cent_ra, cent_dec)
    theta_prior = dict(PDF='exp', max=6., scale=0.5)

    # Pretend numba is absent regardless of the environment
    monkeypatch.setattr(bayesian, 'HAS_NUMBA', False)

    p_ref = bayesian.px_Oi_fixedgrid(
        10., localiz, cand_coords, cand_ang_size, theta_prior,
        use_numba=False)
    with pytest.warns(UserWarning):
        p_fb = bayesian.px_Oi_fixedgrid(
            10., localiz, cand_coords, cand_ang_size, theta_prior,
            use_numba=True)
    assert np.allclose(p_ref, p_fb, rtol=1e-10, atol=1e-15)


@numba_required
def test_px_Oi_fixedgrid_numba_speed_up(capsys):
    """Time numba vs numpy px_Oi_fixedgrid and report the speed-up.

    The first numba call pays JIT-compile time, so the kernel is warmed
    up once before timing.
    """
    cent_ra, cent_dec = 120.0, 32.0
    localiz = _eellipse_localiz(cent_ra, cent_dec, a=1.0, b=0.6, theta=30.)
    cand_coords, cand_ang_size = _make_candidates(
        cent_ra, cent_dec, ncand=50)
    theta_prior = dict(PDF='exp', max=6., scale=0.5)

    def _np(*_):
        return bayesian.px_Oi_fixedgrid(
            10., localiz, cand_coords, cand_ang_size, theta_prior,
            use_numba=False)

    def _nb(*_):
        return bayesian.px_Oi_fixedgrid(
            10., localiz, cand_coords, cand_ang_size, theta_prior,
            use_numba=True)

    # Warm up the JIT (compile) before timing
    p_nb = _nb()
    assert np.allclose(_np(), p_nb, rtol=1e-10, atol=1e-15)

    t_np = _time_call(_np, None, None, None)
    t_nb = _time_call(_nb, None, None, None)

    speedup = t_np / t_nb if t_nb > 0 else np.inf
    with capsys.disabled():
        print(
            "\n  px_Oi_fixedgrid (50 cand, 200x200): "
            "numpy %8.2f ms  numba %8.2f ms  speed-up %5.1fx"
            % (t_np * 1e3, t_nb * 1e3, speedup))
