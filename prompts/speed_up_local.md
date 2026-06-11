# Speeding up the PATH calculations

## Goals

This document is a collection of ideas and prompts for speeding up the PATH calculations for the bayesian.px_Oi_local method.  

## Coding

Here are guidelines for coding: 

- Examine the work described in the speed_up_fixed.md document.
- Use Python where possible
- Add inline comments to explain the effort
- Reuse existing code when possible, especially in bayesian.py and localization.py
- Place import statements at the top of the file.
- Include a description of inputs/outputs in the doc string of all methods
- Use lines of code that are less than 80 characters wide

## Running

If you need to run python, use the "astro14" environment.  

## Testing

As we will be making a number of significant changes to the code, we will need to generate a set of tests to verify the code as we develop the new method.  

1.  Generate a new module named tests_local.py in astropath/tests.  In it place a series of tests that compare the current px_Oi_local method against the px_Oi_fixedgrid method run with a fine grid.  For this:

- Draw on examples in the Notebooks in the astropath/calculations/step_size directory.
- Use a fine enough step size for the px_Oi_fixedgrid method to be considered accurate
- Generate examples for at least these 3 cases:
  * A large galaxy (10") with a small localization (1")
  * A large galaxy (10") with a large localization (10")
  * A small galaxy (1") with a small localization (1")
  * A small galaxy (1") with a large localization (10")
  * A small galaxy (0.5") with an ellipse localization (a=10", b=0.2")
- Log your work in the "Logs" section below.

## Development

We wish to speed up the calculations in bayesian.py:px_Oi_local.  
We will begin by removing all astropy calculations and replacing them with numpy calculations.  We will then use numba to speed up the code.

1. Examine the code in bayesian.py:px_Oi_local.  Please:

- Examine the code
- Identify the bottlenecks
- Figure out how to replace all astropy calculations 
- When we test, we will use the tests_local.py module to verify the code.  It will not be necessary to preserve a copy of the original algorithm
- Log your work in the "Logs" section below.
- Make recommendations 
- Do not generate any new code yet
- If you have any questions, ask them in the Logs and I will clarify at my next prompt

2. Write new code to replace the astropy calculations with numpy calculations. Please:

- Here are answers to your questions:
  * Only consider the ellipse localization type for now.
  * Flat sky is acceptable
- Remember that we intend to use numba to speed up the code.  Consider that in your programming
- Verify using the tests_local.py module that the code is working as expected.
- Add new code to profiling.py to profile this method
- Log your work in the "Logs" section below.

## Profiling

## Docs

## Prompts

1. Read this doc.  Proceed with the 1st item under Testing.
2. Read this doc.  Proceed with the 1st item under Development.
3. Read this doc.  Proceed with the 2nd item under Development.

## Logging

The "Logs" section will record Claude's work.  Please use the following format:

### <Date> (Short summary of the work)

<Detailed description of the work and what you learned>

### <Date> (Short summary of the work)

<Detailed description of the work and what you learned>

...


## Logs

### 2026-06-10 (Generated tests_local.py: px_Oi_local vs fine fixed grid)

**What I built:** A new test module
`astropath/tests/tests_local.py`, modeled after `test_speed_up.py` and
the `calculations/step_size/test_accuracy.py` notebook helpers. It
verifies that `bayesian.px_Oi_local` reproduces
`bayesian.px_Oi_fixedgrid` run on a fine, large-enough grid (taken as
the accurate reference).

**Key insight (why this is a valid comparison):** both methods evaluate
the same integral `p(x|O_i) = ∫ L(w-x) p(w|O_i) dw`, just discretized on
different grids. `px_Oi_local` centers a grid on the *galaxy* sized to
the offset-prior support (`phi*max`); `px_Oi_fixedgrid` centers on the
*transient*. The integrand is non-zero only within `phi*max` of the
galaxy (the offset prior has a hard cutoff there), so a fixed-grid box
of `offset + phi*max + 3"` fully contains it. I confirmed the reference
value is unchanged when the box is enlarged further, so box size is not
limiting.

**Convergence measured:** `px_Oi_local`'s relative error vs the fine
fixed grid is first-order in the (relative) step size, ≈ `0.167 *
step_size`:

| step_size | reldiff |
|-----------|---------|
| 0.100     | -1.67%  |
| 0.050     | -0.83%  |
| 0.020     | -0.33%  |
| 0.010     | -0.17%  |
| 0.005     | -0.08%  |

So the two methods provably converge to the same number.

**Cases covered** (all 5 from the prompt — the prompt says "at least
3" but lists 5):
- large galaxy 10" / small loc 1"
- large galaxy 10" / large loc 10"
- small galaxy 1"  / small loc 1"
- small galaxy 1"  / large loc 10"
- small galaxy 0.5" / ellipse loc a=10", b=0.2"

**Tests written:**
1. `test_local_matches_fine_fixedgrid` (parametrized over the 5 cases) —
   runs the local method at `step_size=0.02` and asserts agreement
   with the fine fixed grid to `rtol=1e-2`. Observed reldiff is a
   remarkably consistent ≈ -3.34e-3 across every regime.
2. `test_local_multi_candidate_matches_fixedgrid` — drives all regimes
   through the candidate loop (a guard for future vectorization/numba),
   each checked against its own reference.
3. `test_local_convergence_to_fixedgrid` — confirms the error decreases
   monotonically and roughly halves when the step halves (first-order).

**Result:** `pytest astropath/tests/tests_local.py -s` → 7 passed in
~10s. The reference grids stay tractable (largest is 3840² for the
10"-galaxy / 1"-loc case, ~3s; the b=0.2" ellipse needs a 0.007" step,
~1s). These now serve as the regression baseline for the upcoming
numpy/numba rewrite of `px_Oi_local`.

**Question for next prompt:** the `Development` section says to remove
astropy from `px_Oi_local`. Should I keep the existing astropy-based
`px_Oi_local` frozen as a reference (as was done for `calc_LWx` /
`px_Oi_fixedgrid` in `test_speed_up.py`), or is the fine-fixed-grid
comparison in `tests_local.py` sufficient as the regression anchor?

> Answer (recorded by user in the next prompt): no need to preserve a
> copy of the original algorithm; `tests_local.py` is the regression
> anchor.

### 2026-06-10 (Analysis of px_Oi_local: bottlenecks + astropy removal)

**Code examined:** `bayesian.px_Oi_local` (lines 347-411).  It loops
over candidates; for each it (1) builds a square grid centered on the
galaxy, sized `box_hwidth = phi*max` with spacing
`step_size_phi = phi*step_size`, (2) computes `theta = sqrt(x^2+y^2)`
and `p_wOi = pw_Oi(theta, phi, prior)`, (3) builds flat-sky `ra/dec`
around the galaxy, (4) calls `localization.calc_LWx(ra, dec, localiz)`,
(5) accumulates `sum(L_wx * p_wOi) * step_size_phi^2`.

**Astropy inventory (everything in the loop body):**
- `for icand, cand_coord in enumerate(cand_coords)` -- iterating a
  SkyCoord array; each item is a scalar SkyCoord (slicing overhead).
- `cand_coord.ra.deg` (line 396).
- `np.cos(cand_coord.dec).value` (line 397) -- Angle trig + `.value`.
- `cand_coord.dec.deg` (line 398).
- Inside `calc_LWx`: `localiz['center_coord'].ra.deg` / `.dec.deg`
  (read once per call -> N times).  The eellipse *math* in calc_LWx is
  already pure numpy (done in the earlier speed_up_fixed.md work).
- No `units` usage in px_Oi_local itself (the module-level
  `from astropy import units` is only used by the deprecated
  `px_Oi_orig`).

**Profile (200 candidates, step_size=0.1, mixed sizes, 3 reps,
cProfile, cumulative):**

| component | share of runtime |
|-----------|------------------|
| `calc_LWx` (numpy Vincenty sep + PA + Gaussian) | ~58% |
| astropy coordinate access/iteration (SkyCoord, Angle, represent_as) | ~25-30% |
| `pw_Oi` | ~5% |
| `np.meshgrid` + grid build | ~3% |

Wall time ~1.0 ms/candidate at step_size=0.1.

**Bottleneck #1 -- `calc_LWx` doing spherical trig on a flat-sky grid.**
px_Oi_local already *builds* `ra/dec` with a flat-sky tangent-plane
approximation, then `calc_LWx` (eellipse) runs the FULL Vincenty
angular separation + position angle (`arctan2`, `hypot`, several
`sin/cos`) to recover offsets it effectively already knows.  Micro-bench
on a 120x120 grid: `calc_LWx` (Vincenty) = 601 us/call vs a direct
flat-sky Gaussian = 100 us/call -- **~6x** -- agreeing to a max
fractional difference of 1.1e-5 (negligible at arcsec grid scale).

**Bottleneck #2 -- per-candidate astropy coordinate access.**  Same
issue fixed earlier in `px_Oi_fixedgrid`: iterating a SkyCoord and
reading `.ra/.dec` per candidate is pure overhead (~25-30% here).

**Key structural finding -- the grid and `p_wOi` are phi-independent.**
Because `box_hwidth = phi*max` and `step_size_phi = phi*step_size`,

    ngrid = round(2*box_hwidth/step_size_phi) = round(2*max/step_size)

so `ngrid` is the SAME for every candidate (e.g. 120 for max=6,
step_size=0.1; verified).  Writing the grid as `phi*max * U` for a
fixed normalized meshgrid `U`, one finds `theta = phi*max*R` and:
- core: `p_wOi = (1/(max*R+1)) / norm(phi)`
- exp:  `p_wOi = exp(-max*R/scale) / norm(phi)`
- uniform: `p_wOi = 1 / norm(phi)`

i.e. the *shape* of `p_wOi` (and the cutoff mask `R<1`) does NOT depend
on phi -- only the scalar `norm(phi)` and the `step_size_phi^2 =
(phi*step_size)^2` factor do.  So `meshgrid`, `theta`, `R`, and the
exp/core grid can be computed ONCE and reused for all candidates.

**Important constraint -- px_Oi_local must stay generic over
localization type.**  Its docstring ("best for large localization
areas which cover a large area of the sky") points at healpix maps
(CHIME baseband), and `calc_LWx` also serves `wcs`.  For healpix/wcs,
`L_wx` is a pixel lookup that genuinely needs `ra/dec`, so we cannot
drop `calc_LWx` wholesale -- the flat-sky shortcut (Bottleneck #1)
applies only to the `eellipse` branch.

**Plan to replace astropy (for the next coding step -- not done yet):**
1. Pre-extract `cand_ra = cand_coords.ra.deg`,
   `cand_dec = cand_coords.dec.deg`,
   `cos_cand_dec = np.cos(np.radians(cand_dec))` ONCE before the loop
   (mirrors `px_Oi_fixedgrid` lines 269-271).  Index these in the loop
   instead of touching the SkyCoord.  This removes ALL astropy from the
   loop body.
2. Pre-extract `center_ra/center_dec` from `localiz['center_coord']`
   once, so the per-call astropy read inside `calc_LWx` is avoided for
   the local method.  (Either compute a flat-sky `L_wx` inline for the
   eellipse type, or pass scalars; see recommendations.)
3. After (1)-(2) the only astropy left is constructing the inputs;
   px_Oi_local's own computation becomes pure numpy.

**Recommendations (ordered by payoff; later items target the numba
step the doc mentions, NOT this one):**
- **R1 (biggest, eellipse):** add a flat-sky `L_wx` path for the
  eellipse type that reuses the rotation/Gaussian already in
  `calc_LWx` but skips Vincenty -- ~6x on the dominant cost.  Keep
  calling `calc_LWx` for healpix/wcs.  Consider factoring the eellipse
  Gaussian into a small numpy helper in `localization.py` so both
  `calc_LWx` and the local fast path share one definition.
- **R2:** precompute the phi-independent normalized grid + `p_wOi`
  shape + cutoff mask once; per candidate only rescale by `norm(phi)`
  and multiply the sum by `(phi*step_size)^2`.  Removes meshgrid/exp
  from the per-candidate loop entirely.
- **R3:** pre-extract candidate/center coordinates (the astropy-removal
  in the plan above).
- **R4 (numba step):** with R1-R3 the per-candidate kernel reduces to
  "offset the fixed grid by the galaxy's center offset, evaluate the
  eellipse Gaussian, multiply by the cached `p_wOi` shape, sum" -- a
  tight numeric loop well suited to an `@njit` kernel analogous to
  `px_Oi_numba`, computing the eellipse `L_wx` inline.

**Questions for the next prompt:**
1. Should the optimization target ONLY the `eellipse` localization
   type (with healpix/wcs falling back to the current `calc_LWx`
   path), or do you also want a fast path for healpix?  (My plan above
   keeps healpix/wcs on `calc_LWx`.)
2. The flat-sky `L_wx` for eellipse differs from the current Vincenty
   result by ~1e-5 fractionally.  Confirm that adopting the flat-sky
   approximation here is acceptable (it is already used to build the
   grid, and `px_Oi_fixedgrid` makes the same approximation).

### 2026-06-10 (Rewrote px_Oi_local in numpy; ~6.5x faster)

**Done** (Development item 2). User answers applied: only the
`eellipse` type gets the fast path; flat-sky is acceptable.

**Changes to `bayesian.px_Oi_local`:**
- **Astropy removed from the loop (R3):** `cand_coords.ra.deg /
  .dec.deg` and `cos(dec)` are pre-extracted to numpy arrays ONCE
  before the loop; the body indexes those arrays instead of iterating
  the SkyCoord.  For the eellipse path the localization center is also
  pulled out once.  No astropy remains in the per-candidate loop.
- **Normalized grid built ONCE (R2, partial):** since `ngrid =
  2*max/step_size` is phi-independent, the loop now reuses a single
  normalized meshgrid (U, V, R) rescaled by `phi*max` per candidate
  instead of rebuilding `np.meshgrid` each time.  (I deliberately did
  NOT do the full analytic phi-cancellation of `norm(phi)` -- it would
  obscure the algorithm and is unnecessary for the numba step, which
  recomputes the PDF per pixel anyway.  `pw_Oi` is still reused as-is.)
- **Flat-sky eellipse L_wx (R1, the big win):** for `eellipse`, L(w-x)
  is now evaluated directly in the tangent plane -- galaxy-center
  offset (E0, N0) plus the grid offsets, rotated into the ellipse frame
  (same `dtheta = 90-PA` convention as `calc_LWx`), then the 2D
  Gaussian.  This skips astropy's Vincenty separation / position angle.
- **Generic fallback preserved:** non-eellipse types (healpix/wcs)
  still build flat-sky `ra/dec` and call `localization.calc_LWx`, so
  production GW/healpix runs are unaffected.  (Caught by the existing
  `test_path.py::test_gw`, which has no `center_coord`; the center
  extraction is now guarded inside the eellipse branch.)
- **numba-readiness:** the eellipse per-candidate body is now a tight
  numeric block (scalars + the fixed normalized grid -> rotation,
  Gaussian, weighted sum) that maps cleanly onto a future `@njit`
  kernel analogous to `px_Oi_numba`.

**Verification:** `pytest astropath/tests/tests_local.py -s` -> 8
passed.  The reldiff vs the fine fixed grid is unchanged from the old
Vincenty version (e.g. -3.345e-3 vs -3.348e-3 for the large-galaxy/
small-loc case) -- i.e. the flat-sky approximation contributes <1e-5,
far below the discretization error and the 1% tolerance.  Full suite
(`test_bayesian`, `test_path`, `test_speed_up`, `tests_local`) passes.

**Speed-up:** in the 200-candidate / step_size=0.1 / eellipse scenario,
0.154 ms/candidate vs ~1.000 ms/candidate before -> **~6.5x** (Vincenty
removal ~6x + per-candidate astropy removal).

**Profiling added** (`astropath/profiling.py`): new
`run_profiling_local()` sweeps `step_size` (per-candidate grid =
2*max/step_size) for the 50-candidate eellipse scenario, plus
`plot_local_results()` (purple curve, x-axis = sqrt(per-candidate
pixels), 10 s reference line) writing `profiling_local_timing.png`.
`main()` now runs and reports both the fixed-grid and local sweeps.
Sample (50 candidates):

| step_size | per-cand grid | total ms | ms/cand |
|-----------|---------------|----------|---------|
| 0.500     | 24x24         | 0.9      | 0.018   |
| 0.250     | 48x48         | 1.8      | 0.035   |
| 0.100     | 120x120       | 7.8      | 0.155   |
| 0.050     | 240x240       | 32.8     | 0.655   |
| 0.025     | 480x480       | 182.9    | 3.659   |
