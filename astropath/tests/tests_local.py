""" Accuracy tests for bayesian.px_Oi_local.

These tests establish a correctness baseline for ``px_Oi_local`` BEFORE
we start optimizing it (numpy rewrite, then numba; see
``prompts/speed_up_local.md``).  Both ``px_Oi_local`` and
``px_Oi_fixedgrid`` evaluate the same integral

    p(x|O_i) = integral over w of  L(w-x) * p(w|O_i) dw

just discretized on different grids: ``px_Oi_local`` builds a grid
centered on the galaxy and sized to the offset-prior support
(phi*max), while ``px_Oi_fixedgrid`` builds a grid centered on the
transient.  Run with a sufficiently fine step the fixed-grid result is
effectively the truth, so we use it as the reference and check that the
(coarser) local method reproduces it.

Empirically the local method's RELATIVE error is ~0.17*step_size (first
order in the grid spacing): step_size=0.1 -> ~1.7%, 0.05 -> ~0.8%,
0.02 -> ~0.33%.  The tests below run the local method at step_size=0.02
and require agreement to better than 1% with the fine fixed grid.

The cases requested in the prompt cover the relevant size regimes:
  * large galaxy (10") with a small localization (1")
  * large galaxy (10") with a large localization (10")
  * small galaxy (1")  with a small localization (1")
  * small galaxy (1")  with a large localization (10")
  * small galaxy (0.5") with an ellipse localization (a=10", b=0.2")
"""

import numpy as np

from astropy.coordinates import SkyCoord
from astropy import units

from astropath import bayesian

import pytest

# Offset prior shared by every comparison.  Exponential PDF with a hard
# cutoff at max*phi (matches the priors used elsewhere in the package).
THETA_PRIOR = dict(max=6., PDF='exp', scale=0.5)

# Reference transient = localization center (mid-declination ICRS so the
# cos(dec) flat-sky scaling is exercised but not degenerate).
CENT_RA, CENT_DEC = 120.0, 32.0

# Local-method grid step used for the matching tests.  At this (relative)
# step the local error is ~0.33%, comfortably inside RTOL below.
LOCAL_STEP = 0.02
RTOL = 1.0e-2

# Comparison cases.  Each tuple is:
#   (label, phi_gal, a, b, pa_ee, gal_offset, gal_pa)
#     phi_gal    : galaxy angular size (arcsec)
#     a, b       : localization ellipse semi-axes (arcsec); a==b -> circle
#     pa_ee      : position angle of the localization ellipse (deg)
#     gal_offset : galaxy offset from the transient center (arcsec)
#     gal_pa     : position angle of that offset (deg E of N)
CASES = [
    ("large_gal10_small_loc1",  10.,  1.,  1.,  0., 1.0,  0.),
    ("large_gal10_large_loc10", 10., 10., 10.,  0., 3.0,  0.),
    ("small_gal1_small_loc1",    1.,  1.,  1.,  0., 1.0,  0.),
    ("small_gal1_large_loc10",   1., 10., 10.,  0., 3.0,  0.),
    ("small_gal0p5_ellipse",     0.5, 10., 0.2, 0., 2.0, 90.),
]


def _center():
    """Return the transient/localization center coordinate.

    Returns:
        SkyCoord: ICRS center coordinate (CENT_RA, CENT_DEC).
    """
    return SkyCoord(ra=CENT_RA, dec=CENT_DEC, unit='deg')


def _eellipse_localiz(a, b, pa_ee):
    """Build an eellipse localization dict centered on the transient.

    Args:
        a (float): Semi-major axis of the error ellipse, arcsec.
        b (float): Semi-minor axis of the error ellipse, arcsec.
        pa_ee (float): Position angle of the ellipse, deg (E of N).

    Returns:
        dict: Localization dict (eellipse type) with center_coord set.
    """
    return dict(type='eellipse',
                center_coord=_center(),
                eellipse=dict(a=a, b=b, theta=pa_ee))


def _galaxy(offset, gal_pa):
    """Build a single-candidate SkyCoord offset from the transient.

    Args:
        offset (float): Angular offset from the center, arcsec.
        gal_pa (float): Position angle of the offset, deg (E of N).

    Returns:
        SkyCoord: Length-1 SkyCoord array for the candidate galaxy.
    """
    gal = _center().directional_offset_by(
        gal_pa * units.deg, offset * units.arcsec)
    # Re-wrap as a length-1 array so it iterates like the real input.
    return SkyCoord([gal.ra.deg], [gal.dec.deg], unit='deg')


def _fine_fixedgrid_reference(phi, a, b, pa_ee, offset, gal_pa):
    """Accurate p(x|O_i) from px_Oi_fixedgrid on a fine, large-enough grid.

    The integrand L(w-x)*p(w|O_i) is non-zero only within phi*max of the
    galaxy (the offset prior's hard cutoff), so a box reaching
    ``offset + phi*max`` from the transient (plus a small margin) fully
    contains it.  The step resolves the smallest relevant scale
    (the smaller of the ellipse minor axis and the galaxy size), making
    the fixed-grid result effectively exact for our purposes.

    Args:
        phi (float): Galaxy angular size, arcsec.
        a (float): Localization semi-major axis, arcsec.
        b (float): Localization semi-minor axis, arcsec.
        pa_ee (float): Localization ellipse PA, deg.
        offset (float): Galaxy offset from the transient, arcsec.
        gal_pa (float): PA of the galaxy offset, deg.

    Returns:
        float: Reference p(x|O_i) for the single candidate.
    """
    localiz = _eellipse_localiz(a, b, pa_ee)
    cand = _galaxy(offset, gal_pa)
    cand_ang_size = np.array([phi])
    # Box must cover the offset-prior support around the galaxy.
    box_hwidth = offset + phi * THETA_PRIOR['max'] + 3.0  # arcsec
    # Resolve the smallest scale; floor keeps grids tractable.
    fine_step = max(min(b, phi) / 30., 0.005)             # arcsec
    p_ref = bayesian.px_Oi_fixedgrid(
        box_hwidth, localiz, cand, cand_ang_size,
        THETA_PRIOR, step_size=fine_step)
    return p_ref[0]


@pytest.mark.parametrize(
    "label,phi,a,b,pa_ee,offset,gal_pa", CASES,
    ids=[c[0] for c in CASES])
def test_local_matches_fine_fixedgrid(
        label, phi, a, b, pa_ee, offset, gal_pa, capsys):
    """px_Oi_local must match the fine fixed-grid reference (<1%).

    Runs the local method at step_size=LOCAL_STEP and compares to the
    accurate fixed-grid value for each size regime.
    """
    localiz = _eellipse_localiz(a, b, pa_ee)
    cand = _galaxy(offset, gal_pa)
    cand_ang_size = np.array([phi])

    p_ref = _fine_fixedgrid_reference(
        phi, a, b, pa_ee, offset, gal_pa)
    p_loc = bayesian.px_Oi_local(
        localiz, cand, cand_ang_size, THETA_PRIOR,
        step_size=LOCAL_STEP)[0]

    rel = (p_loc - p_ref) / p_ref
    with capsys.disabled():
        print("\n  %-24s local=%11.5e  fixed(fine)=%11.5e  "
              "reldiff=%+.3e" % (label, p_loc, p_ref, rel))

    assert np.isfinite(p_loc) and p_loc > 0
    assert np.isclose(p_loc, p_ref, rtol=RTOL, atol=0.0)


# Cases where the localization minor axis is smaller than the galaxy
# (b < phi) -- these trigger the localization-centered L_wx correction.
CORR_CASES = [c for c in CASES if c[3] < c[1]]  # b (idx 3) < phi (idx 1)


@pytest.mark.parametrize(
    "label,phi,a,b,pa_ee,offset,gal_pa", CORR_CASES,
    ids=[c[0] for c in CORR_CASES])
@pytest.mark.parametrize("coarse_step", [0.05, 0.1])
def test_local_correction_coarse_step(
        label, phi, a, b, pa_ee, offset, gal_pa, coarse_step, capsys):
    """When b < phi, the _Lwx_correction keeps px_Oi_local accurate even
    at coarse (galaxy-relative) step sizes.

    The fast aligned-grid correction divides the under-resolved raw sum
    by the discrete "total L_wx"; the aliasing cancels, recovering ~0.1-
    0.2% accuracy where the UNcorrected raw sum is biased by ~1-2% (the
    O(step) error).  Tolerance is well below that raw bias but looser
    than the old fine-grid helper (the new method is ~1%-class by
    design).
    """
    localiz = _eellipse_localiz(a, b, pa_ee)
    cand = _galaxy(offset, gal_pa)
    cand_ang_size = np.array([phi])

    p_ref = _fine_fixedgrid_reference(
        phi, a, b, pa_ee, offset, gal_pa)
    p_loc = bayesian.px_Oi_local(
        localiz, cand, cand_ang_size, THETA_PRIOR,
        step_size=coarse_step)[0]

    rel = (p_loc - p_ref) / p_ref
    with capsys.disabled():
        print("\n  %-24s step=%.2f local=%11.5e fixed(fine)=%11.5e "
              "reldiff=%+.3e" % (label, coarse_step, p_loc, p_ref, rel))

    # Much tighter than the generic RTOL: the correction removes the
    # bulk of the step-size bias for these cases.
    assert np.isclose(p_loc, p_ref, rtol=5.0e-3, atol=0.0)


# Very small (sub-arcsec) circular localizations -- the deeply
# under-resolved regime (grid spacing phi*step >> b).  phi is kept <= 1"
# so the fine fixed-grid reference stays well under 5000 cells per side.
SMALL_LOC_CASES = [
    ("small_loc_gal0p3", 0.3, 0.1, 0.1, 0., 0.3, 0.),
    ("small_loc_gal0p6", 0.6, 0.1, 0.1, 0., 0.3, 0.),
    ("small_loc_gal1",   1.0, 0.1, 0.1, 0., 0.3, 0.),
]


@pytest.mark.parametrize(
    "label,phi,a,b,pa_ee,offset,gal_pa", SMALL_LOC_CASES,
    ids=[c[0] for c in SMALL_LOC_CASES])
def test_local_small_localization(
        label, phi, a, b, pa_ee, offset, gal_pa, capsys):
    """px_Oi_local stays accurate for a very small localization (0.1").

    Here ``b < phi`` so the _Lwx_correction fires, and the galaxy grid
    badly under-resolves the 0.1" localization (spacing phi*step can be
    many times b).  The correction must still recover the fine fixed-grid
    value.  Run at the default step (0.05).
    """
    localiz = _eellipse_localiz(a, b, pa_ee)
    cand = _galaxy(offset, gal_pa)
    cand_ang_size = np.array([phi])

    p_ref = _fine_fixedgrid_reference(
        phi, a, b, pa_ee, offset, gal_pa)
    p_loc = bayesian.px_Oi_local(
        localiz, cand, cand_ang_size, THETA_PRIOR, step_size=0.05)[0]

    rel = (p_loc - p_ref) / p_ref
    with capsys.disabled():
        print("\n  %-18s local=%11.5e  fixed(fine)=%11.5e  "
              "reldiff=%+.3e" % (label, p_loc, p_ref, rel))

    assert np.isfinite(p_loc) and p_loc > 0
    assert np.isclose(p_loc, p_ref, rtol=RTOL, atol=0.0)


def test_local_multi_candidate_matches_fixedgrid(capsys):
    """All cases at once: px_Oi_local on a multi-candidate input.

    Drives the candidate loop with every regime in a single call (a
    useful guard for future vectorization/numba work) and checks each
    against its fine fixed-grid reference.
    """
    # Build one SkyCoord holding all candidate galaxies + their sizes.
    ras, decs, sizes, refs = [], [], [], []
    for label, phi, a, b, pa_ee, offset, gal_pa in CASES:
        gal = _galaxy(offset, gal_pa)
        ras.append(gal.ra.deg[0])
        decs.append(gal.dec.deg[0])
        sizes.append(phi)
        refs.append(_fine_fixedgrid_reference(
            phi, a, b, pa_ee, offset, gal_pa))
    cand = SkyCoord(ra=ras, dec=decs, unit='deg')
    cand_ang_size = np.array(sizes)
    refs = np.array(refs)

    # The localization differs per case, so px_Oi_local cannot be run
    # once for all candidates; we drive it per candidate but through the
    # same array-style call signature used in production.
    p_loc = np.array([
        bayesian.px_Oi_local(
            _eellipse_localiz(a, b, pa_ee),
            SkyCoord([ra], [dec], unit='deg'),
            np.array([phi]), THETA_PRIOR, step_size=LOCAL_STEP)[0]
        for (label, phi, a, b, pa_ee, offset, gal_pa), ra, dec
        in zip(CASES, ras, decs)])

    rel = (p_loc - refs) / refs
    with capsys.disabled():
        print("\n  multi-candidate max |reldiff| = %.3e" %
              np.max(np.abs(rel)))

    assert np.all(np.isfinite(p_loc)) and np.all(p_loc > 0)
    assert np.allclose(p_loc, refs, rtol=RTOL, atol=0.0)


def test_local_convergence_to_fixedgrid():
    """The local error must shrink ~linearly as step_size decreases.

    Uses the small-galaxy/small-localization case (cheap, well behaved)
    and confirms the relative error vs the fine fixed grid both (a)
    decreases monotonically and (b) roughly halves when the step is
    halved -- i.e. first-order convergence to the same integral.
    """
    phi, a, b, pa_ee, offset, gal_pa = 1., 1., 1., 0., 1.0, 0.
    localiz = _eellipse_localiz(a, b, pa_ee)
    cand = _galaxy(offset, gal_pa)
    cand_ang_size = np.array([phi])

    p_ref = _fine_fixedgrid_reference(
        phi, a, b, pa_ee, offset, gal_pa)

    steps = [0.1, 0.05, 0.025]
    rel_err = []
    for ss in steps:
        p_loc = bayesian.px_Oi_local(
            localiz, cand, cand_ang_size, THETA_PRIOR,
            step_size=ss)[0]
        rel_err.append(abs(p_loc - p_ref) / p_ref)
    rel_err = np.array(rel_err)

    # Monotonic decrease with finer steps.
    assert np.all(np.diff(rel_err) < 0)
    # First-order: halving the step should roughly halve the error.
    # Allow a generous band around the ideal factor of 2.
    for i in range(len(steps) - 1):
        ratio = rel_err[i] / rel_err[i + 1]
        assert 1.6 < ratio < 2.4
