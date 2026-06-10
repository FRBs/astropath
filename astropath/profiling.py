""" Profiling utilities for the PATH calculations.

This module times the two core PATH steps that were targeted for the
numpy speed-up work (see prompts/speed_up.md):

  * ``localization.calc_LWx`` -- the localization term L(w-x)
  * ``bayesian.px_Oi_fixedgrid`` -- the fixed-grid p(x|O_i)

It is modeled on the ``calculations/step_size/Profiling.ipynb``
notebook: a faux FRB, a circular error ellipse, and a couple of
candidate galaxies.  Rather than timing each sub-step interactively, it
sweeps a range of grid sizes (by varying ``step_size``), builds a table
of timings, prints it, and writes a figure.

Run from the command line::

    python -m astropath.profiling

or import and call :func:`run_profiling` directly.
"""
import os
import time

import numpy as np
import pandas

import matplotlib
matplotlib.use('Agg')  # headless-safe; no interactive display needed
import matplotlib.pyplot as plt  # noqa: E402

from astropy.coordinates import SkyCoord
from astropy import units

from astropath import localization
from astropath import bayesian


# Faux FRB and analysis defaults (mirroring the Profiling notebook)
FRB_RADEC = '21h44m25.255s -40d54m00.10s'
BOX_HWIDTH = 90.  # arcsec; half-width of the analysis box
# step_size sweep (arcsec).  ngrid = 2*box_hwidth/step_size, so with
# box_hwidth=90 these give side lengths of 360 ... 7200 pixels.
DEFAULT_STEP_SIZES = [0.5, 0.25, 0.1, 0.05, 0.025]


NCAND = 50  # number of candidate galaxies in the profiling scenario


def default_setup(ncand=NCAND):
    """Build the faux-FRB profiling scenario.

    Mirrors calculations/step_size/Profiling.ipynb (circular 5" error
    ellipse), but scatters ``ncand`` candidate galaxies over a range of
    locations and angular sizes so the per-candidate loop in
    ``px_Oi_fixedgrid`` is exercised realistically.

    Args:
        ncand (int, optional): Number of candidate galaxies.

    Returns:
        tuple: (localiz, cand_coords, cand_ang_size, theta_prior)
            localiz (dict): eellipse localization dict.
            cand_coords (SkyCoord): candidate host coordinates.
            cand_ang_size (np.ndarray): candidate angular sizes, arcsec.
            theta_prior (dict): offset-prior parameters.
    """
    frb_coord = SkyCoord(FRB_RADEC, frame='icrs')
    # Small but not tiny, circular error ellipse
    eellipse = dict(a=5., b=5., theta=0.)
    localiz = dict(type='eellipse',
                   center_coord=frb_coord,
                   eellipse=eellipse)
    # Scatter the candidates deterministically (no RNG) over a range of
    # position angles and separations from the transient.  Separations
    # span 0.5"-10" so some fall well inside the localization and some
    # near its edge.
    idx = np.arange(ncand)
    pa = (idx * 360. / ncand) * units.deg          # spread in PA
    sep = np.linspace(0.5, 10., ncand) * units.arcsec  # spread in offset
    cand_coords = frb_coord.directional_offset_by(pa, sep)
    # Range of angular sizes, 0.2"-3.0"
    cand_ang_size = np.linspace(0.2, 3.0, ncand)  # arcsec
    theta_prior = dict(max=6., PDF='exp', scale=1.)
    return localiz, cand_coords, cand_ang_size, theta_prior


def _time_call(func, args, npix, reps=None):
    """Time a single callable, returning the best wall time in seconds.

    The number of repetitions is scaled down for large grids to keep
    the overall sweep fast, unless an explicit ``reps`` is given.

    Args:
        func (callable): Function to time.
        args (tuple): Positional arguments passed to ``func``.
        npix (int): Number of grid pixels (used to pick repetitions).
        reps (int, optional): Force a fixed number of repetitions.
            Use reps=1 to time a single call (e.g. to *include* the
            numba JIT compilation cost in the first measurement).

    Returns:
        float: Best-of-reps wall-clock time, seconds.
    """
    if reps is None:
        # Fewer reps for big grids (they are slow and stable)
        if npix > 1_000_000:
            reps = 1
        elif npix > 100_000:
            reps = 2
        else:
            reps = 3
    best = np.inf
    for _ in range(reps):
        t0 = time.perf_counter()
        func(*args)
        best = min(best, time.perf_counter() - t0)
    return best


def _build_grid(localiz, box_hwidth, step_size):
    """Build the fixed ra/dec grid the same way px_Oi_fixedgrid does.

    Args:
        localiz (dict): eellipse localization dict (needs center_coord).
        box_hwidth (float): Half-width of the analysis box, arcsec.
        step_size (float): Grid step size, arcsec.

    Returns:
        tuple: (ra, dec, ngrid) with ra/dec numpy arrays in deg and
            ngrid the number of pixels per side.
    """
    ngrid = int(np.round(2 * box_hwidth / step_size))
    x = np.linspace(-box_hwidth, box_hwidth, ngrid)
    xcoord, ycoord = np.meshgrid(x, x)
    # Flat-sky; RA increases in x
    center = localiz['center_coord']
    cos_dec = np.cos(np.radians(center.dec.deg))
    ra = center.ra.deg + xcoord / 3600. / cos_dec
    dec = center.dec.deg + ycoord / 3600.
    return ra, dec, ngrid


def run_profiling(step_sizes=None, box_hwidth=BOX_HWIDTH):
    """Profile calc_LWx and px_Oi_fixedgrid over a range of grid sizes.

    Args:
        step_sizes (list, optional): Grid step sizes (arcsec) to sweep.
            Defaults to :data:`DEFAULT_STEP_SIZES`.
        box_hwidth (float, optional): Analysis-box half-width, arcsec.

    Returns:
        pandas.DataFrame: One row per step size with columns
            ``step_size``, ``ngrid``, ``n_pixels``, ``calc_LWx_s``,
            ``px_Oi_fixedgrid_s`` (numpy) and ``px_Oi_numba_s`` (the
            numba ``use_numba=True`` path; NaN if numba is unavailable).
    """
    if step_sizes is None:
        step_sizes = DEFAULT_STEP_SIZES
    localiz, cand_coords, cand_ang_size, theta_prior = default_setup()

    # Default correction applied to all px_Oi_fixedgrid timings here
    correction = 'p_wO'

    print("Starting profiling sweep over %d step sizes, %d candidates"
          % (len(step_sizes), len(cand_ang_size)))
    print("  using correction=%r" % correction)
    # NOTE: the numba JIT is intentionally NOT pre-warmed.  Its
    # compilation cost is included in the first (smallest-grid) numba
    # timing below, so the reported numba time reflects a realistic
    # one-off sandbox call.

    rows = []
    for step_size in step_sizes:
        # Grid for the standalone calc_LWx timing
        ra, dec, ngrid = _build_grid(localiz, box_hwidth, step_size)
        npix = ngrid * ngrid
        print("  step_size=%.4f  grid=%dx%d (%d pix) ..."
              % (step_size, ngrid, ngrid, npix))

        # Time calc_LWx (localization term only)
        t_lwx = _time_call(localization.calc_LWx,
                           (ra, dec, localiz), npix)
        print("    calc_LWx: %.1f ms" % (t_lwx * 1e3))

        # Time the full fixed-grid p(x|O_i) -- numpy path
        t_pxoi = _time_call(
            lambda *a: bayesian.px_Oi_fixedgrid(
                box_hwidth, localiz, cand_coords, cand_ang_size,
                theta_prior, step_size=step_size, correction=correction),
            (), npix)
        print("    px_Oi_fixedgrid (numpy): %.1f ms" % (t_pxoi * 1e3))

        # Time the numba path (use_numba=True), if available.  reps=1 so
        # the first call's JIT compilation is counted in the timing.
        if bayesian.HAS_NUMBA:
            t_numba = _time_call(
                lambda *a: bayesian.px_Oi_fixedgrid(
                    box_hwidth, localiz, cand_coords, cand_ang_size,
                    theta_prior, step_size=step_size, use_numba=True,
                    correction=correction),
                (), npix, reps=1)
            print("    px_Oi_fixedgrid (numba): %.1f ms (%.1fx)"
                  % (t_numba * 1e3, t_pxoi / t_numba))
        else:
            t_numba = np.nan

        rows.append(dict(step_size=step_size, ngrid=ngrid,
                         n_pixels=npix, calc_LWx_s=t_lwx,
                         px_Oi_fixedgrid_s=t_pxoi,
                         px_Oi_numba_s=t_numba))

    print("Profiling sweep complete.")
    return pandas.DataFrame(rows)


def plot_results(df, outfile):
    """Plot the timing results vs grid size on log-log axes.

    Args:
        df (pandas.DataFrame): Output of :func:`run_profiling`.
        outfile (str): Path to write the PNG figure.

    Returns:
        str: The path the figure was written to.
    """
    # x-axis in sqrt(pixels) = grid side length
    sqrt_pix = np.sqrt(df['n_pixels'])
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(sqrt_pix, df['calc_LWx_s'], 'o-',
            label='calc_LWx')
    ax.plot(sqrt_pix, df['px_Oi_fixedgrid_s'], 's-', color='red',
            label='px_Oi_fixedgrid (numpy)')
    # numba curve (only if it was timed)
    if 'px_Oi_numba_s' in df and df['px_Oi_numba_s'].notna().any():
        ax.plot(sqrt_pix, df['px_Oi_numba_s'], '^-', color='green',
                label='px_Oi_fixedgrid (numba)')
    # Reference line at 10 s
    ax.axhline(10., color='dimgray', linestyle='--', linewidth=1.5)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('sqrt(grid pixels)')
    ax.set_ylabel('Time (s)')
    ax.set_title('PATH profiling')
    ax.grid(True, which='both', alpha=0.3)
    ax.set_xlim(None, 10000.)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outfile, dpi=120)
    plt.close(fig)
    return outfile


def main():
    """Run the profiling sweep, print the table, and save the figure.

    Args:
        None

    Returns:
        pandas.DataFrame: The timing table (also printed to screen).
    """
    df = run_profiling()

    # Print a readable table (times in ms for legibility)
    show = df.copy()
    show['calc_LWx_ms'] = show['calc_LWx_s'] * 1e3
    show['px_Oi_numpy_ms'] = show['px_Oi_fixedgrid_s'] * 1e3
    show['px_Oi_numba_ms'] = show['px_Oi_numba_s'] * 1e3
    # numba speed-up factor over the numpy path
    show['numba_speedup'] = (show['px_Oi_fixedgrid_s']
                             / show['px_Oi_numba_s'])
    cols = ['step_size', 'ngrid', 'n_pixels', 'calc_LWx_ms',
            'px_Oi_numpy_ms', 'px_Oi_numba_ms', 'numba_speedup']
    print('\nPATH profiling results (best-of-reps):')
    print(show[cols].to_string(index=False,
                               float_format=lambda v: '%.3f' % v))

    # Save the figure next to this module
    outfile = os.path.join(os.path.dirname(__file__),
                           'profiling_timing.png')
    plot_results(df, outfile)
    print('\nFigure written to: %s' % outfile)
    return df


if __name__ == '__main__':
    main()
