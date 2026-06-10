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


def default_setup():
    """Build the faux-FRB profiling scenario.

    Mirrors calculations/step_size/Profiling.ipynb: a circular 5"
    error ellipse and two candidate galaxies 1" from the transient.

    Args:
        None

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
    # Two candidates: 1" North and 1" East of the transient
    gal1 = frb_coord.directional_offset_by(0. * units.deg,
                                           1.0 * units.arcsec)
    gal2 = frb_coord.directional_offset_by(90. * units.deg,
                                           1.0 * units.arcsec)
    cand_coords = SkyCoord([gal1.ra, gal2.ra], [gal1.dec, gal2.dec])
    cand_ang_size = np.array([1.0, 0.2])  # arcsec
    theta_prior = dict(max=6., PDF='exp', scale=1.)
    return localiz, cand_coords, cand_ang_size, theta_prior


def _time_call(func, args, npix):
    """Time a single callable, returning the best wall time in seconds.

    The number of repetitions is scaled down for large grids to keep
    the overall sweep fast.

    Args:
        func (callable): Function to time.
        args (tuple): Positional arguments passed to ``func``.
        npix (int): Number of grid pixels (used to pick repetitions).

    Returns:
        float: Best-of-reps wall-clock time, seconds.
    """
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
            ``step_size``, ``ngrid``, ``n_pixels``, ``calc_LWx_s`` and
            ``px_Oi_fixedgrid_s`` (times in seconds).
    """
    if step_sizes is None:
        step_sizes = DEFAULT_STEP_SIZES
    localiz, cand_coords, cand_ang_size, theta_prior = default_setup()

    rows = []
    for step_size in step_sizes:
        # Grid for the standalone calc_LWx timing
        ra, dec, ngrid = _build_grid(localiz, box_hwidth, step_size)
        npix = ngrid * ngrid

        # Time calc_LWx (localization term only)
        t_lwx = _time_call(localization.calc_LWx,
                           (ra, dec, localiz), npix)

        # Time the full fixed-grid p(x|O_i)
        t_pxoi = _time_call(
            bayesian.px_Oi_fixedgrid,
            (box_hwidth, localiz, cand_coords, cand_ang_size,
             theta_prior, step_size), npix)

        rows.append(dict(step_size=step_size, ngrid=ngrid,
                         n_pixels=npix, calc_LWx_s=t_lwx,
                         px_Oi_fixedgrid_s=t_pxoi))

    return pandas.DataFrame(rows)


def plot_results(df, outfile):
    """Plot the timing results vs grid size on log-log axes.

    Args:
        df (pandas.DataFrame): Output of :func:`run_profiling`.
        outfile (str): Path to write the PNG figure.

    Returns:
        str: The path the figure was written to.
    """
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(df['n_pixels'], df['calc_LWx_s'], 'o-',
            label='calc_LWx')
    ax.plot(df['n_pixels'], df['px_Oi_fixedgrid_s'], 's-',
            label='px_Oi_fixedgrid')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Grid size (pixels)')
    ax.set_ylabel('Time (s)')
    ax.set_title('PATH profiling: numpy implementations')
    ax.grid(True, which='both', alpha=0.3)
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
    show['px_Oi_fixedgrid_ms'] = show['px_Oi_fixedgrid_s'] * 1e3
    cols = ['step_size', 'ngrid', 'n_pixels',
            'calc_LWx_ms', 'px_Oi_fixedgrid_ms']
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
