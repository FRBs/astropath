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

- ``bayesian.px_Oi_local`` (the local-grid :math:`p(x|O_i)`, which
  builds one grid per candidate) has likewise been rewritten in pure
  numpy for the ``eellipse`` localization: it pre-extracts the
  coordinates once and evaluates :math:`L(w-x)` directly in a flat-sky
  tangent plane, removing every ``astropy`` call from the per-candidate
  loop. Other localization types (``healpix``, ``wcs``) fall back to
  ``localization.calc_LWx`` as before. See `Local-grid posterior`_
  below for details.

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

Local-grid posterior
====================

``bayesian.px_Oi_local`` evaluates :math:`p(x|O_i)` on a separate grid
per candidate, centered on the galaxy and sized to the offset prior
(``box_hwidth = phi*max``). It is the method of choice for
localizations that span a large area of sky, where
a single fixed grid would be prohibitively large.

For the ``eellipse`` localization the calculation is pure numpy and
flat-sky:

- The per-candidate grid is built once in normalized units and merely
  rescaled per galaxy — the pixel count ``ngrid = 2*max/step_size`` is
  the same for every candidate.

- :math:`L(w-x)` is evaluated directly in the tangent plane (offsets
  rotated into the ellipse frame), matching the spherical ``calc_LWx``
  to ~1e-5 fractionally at arcsec scales.

Here ``step_size`` is *relative* to the galaxy size (default ``0.05``),
so the grid spacing is ``phi*step_size``.

.. note::

   ``px_Oi_local`` is **pure numpy and does not require (or use)
   numba** — for now. The optional numba acceleration described above
   applies only to ``px_Oi_fixedgrid``. ``px_Oi_local`` is fast on its
   own because each per-candidate grid is small.

Small-localization correction
------------------------------

When the localization minor axis :math:`b` is smaller than the galaxy
angular size :math:`\phi`, the galaxy-centered grid under-resolves the
sharp localization and the raw sum is biased low. In that case
``px_Oi_local`` divides the result by a correction factor computed by
``bayesian._Lwx_correction``: the discrete "total :math:`L(w-x)`" on a
small grid that is centered on the localization and *aligned to the
galaxy grid* (same spacing, shifted by an integer number of cells).
Because the localization is sampled at the same sub-cell phase in the
raw sum and in this factor, the under-resolution bias cancels in the
ratio (accurate to ~1%). This is the local analogue of
``px_Oi_fixedgrid``'s ``correction='L_wx'``. The correction grid is
bounded (it is skipped when it would exceed ~5000 cells per side, which
only happens when the localization is already well resolved and no
correction is needed), so it never allocates a large array.

Profiling module
================

The ``astropath.profiling`` module measures these calculations across a
range of grid sizes. It mirrors the
``calculations/step_size/Profiling.ipynb`` notebook: a faux FRB, a
circular error ellipse, and a set of candidate galaxies spread over a
range of locations and angular sizes.

Run it from the command line::

    python -m astropath.profiling

This sweeps a range of ``step_size`` values (hence grid sizes) and
profiles both posterior methods:

- ``calc_LWx``, the numpy ``px_Oi_fixedgrid``, and (if ``numba`` is
  installed) the numba path with its speed-up factor — written to
  ``profiling_timing.png``;

- ``px_Oi_local`` (via ``run_profiling_local``), whose per-candidate
  grid size is ``2*max/step_size`` — written to
  ``profiling_local_timing.png``.

Both timing tables are printed to the screen.

You can also import and call the pieces directly::

    from astropath import profiling

    df = profiling.run_profiling()            # fixed-grid (+ numba)
    profiling.plot_results(df, 'timing.png')

    df_local = profiling.run_profiling_local()   # local-grid method
    profiling.plot_local_results(df_local, 'local.png')

API Reference
=============

.. automodule:: astropath.profiling
   :members:
   :undoc-members:
   :show-inheritance:
