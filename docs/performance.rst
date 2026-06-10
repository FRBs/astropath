***********
Performance
***********

This doc describes the performance-oriented pieces of *astropath*:
the numpy-accelerated core calculations, the optional ``numba``
acceleration for the fixed-grid posterior, and the ``profiling``
module used to measure them.

numpy-accelerated core
======================

The two most expensive steps in a PATH analysis have been
re-implemented in pure numpy:

- ``localization.calc_LWx`` (the localization term :math:`L(w-x)`)
  for the error-ellipse (``eellipse``) localization now computes the
  angular separation and position angle with closed-form numpy
  formulae instead of ``astropy`` ``SkyCoord`` operations. Results are
  numerically identical (to roundoff) and the routine is several times
  faster on large grids.

- ``bayesian.px_Oi_fixedgrid`` (the fixed-grid :math:`p(x|O_i)`)
  pre-extracts the candidate and center coordinates into numpy arrays
  once, rather than reading ``astropy`` attributes inside the
  per-candidate loop.

These changes are transparent: you do not need to do anything to
benefit from them, and the public APIs are unchanged.

Optional numba acceleration
===========================

The per-candidate loop of ``bayesian.px_Oi_fixedgrid`` can optionally
be evaluated with a `numba <https://numba.pydata.org/>`_ kernel
(``bayesian.px_Oi_numba``) that fuses the offset, offset-PDF, product
with :math:`L(w-x)`, and grid sum into a single pass with no full-grid
temporaries.

To enable it, pass ``use_numba=True``::

    from astropath import bayesian

    p_xOi = bayesian.px_Oi_fixedgrid(
        box_hwidth, localiz, cand_coords, cand_ang_size,
        theta_prior, step_size=0.1, use_numba=True)

Key points:

- **Optional.** ``numba`` need *not* be installed. If it is absent (or
  if you simply leave ``use_numba=False``, the default), the
  calculation runs via the standard numpy path. With ``use_numba=True``
  but ``numba`` not installed, the code emits a warning and falls back
  to numpy — it never errors.

- **Default is off.** ``use_numba`` defaults to ``False``, so existing
  code and results are unchanged unless you opt in.

- **Only ``px_Oi_fixedgrid``.** The numba option is available *only*
  for the fixed-grid method. ``px_Oi_local`` and the other routines are
  unaffected. The fused kernel returns only the scalar posterior, so it
  is bypassed (numpy path used) when ``return_grids`` or
  ``return_debug`` is requested.

- **Primarily for sandbox analyses.** numba is recommended mainly for
  interactive/sandbox work and large-grid experiments, where its
  speed-up is largest. The first call pays a one-time JIT compilation
  cost; the win grows with the grid size and the number of candidates
  (roughly 1.5x on small grids up to ~5-6x on a 7200x7200 grid with
  many candidates).

.. note::

   The numba flag is exposed on ``bayesian.px_Oi_fixedgrid`` directly.
   The high-level ``path.PATH.calc_posteriors`` does not currently
   forward ``use_numba``; call ``bayesian.px_Oi_fixedgrid`` directly to
   use it.

Profiling module
================

The ``astropath.profiling`` module measures these calculations across a
range of grid sizes. It mirrors the
``calculations/step_size/Profiling.ipynb`` notebook: a faux FRB, a
circular error ellipse, and a set of candidate galaxies spread over a
range of locations and angular sizes.

Run it from the command line::

    python -m astropath.profiling

This sweeps a range of ``step_size`` values (hence grid sizes), prints
a timing table comparing ``calc_LWx``, the numpy ``px_Oi_fixedgrid``,
and (if ``numba`` is installed) the numba path with its speed-up
factor, and writes a figure (``profiling_timing.png``) of the timings
versus grid size.

You can also import and call it directly::

    from astropath import profiling

    df = profiling.run_profiling()      # returns a pandas DataFrame
    profiling.plot_results(df, 'timing.png')

API Reference
=============

.. automodule:: astropath.profiling
   :members:
   :undoc-members:
   :show-inheritance:
