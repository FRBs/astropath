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

3. I have reduced the default step size to 0.05, relative to the galaxy angular size.  This will not be small enough when the localizatoin error is very small.  In that case, we need to calculate a correction similar to what is done in the fixed method.  The subtlety is that we need to correct it centered on the localization and not at the galaxy.  Please update the code to implement this.  It will mean:

- Checking to see if b is less than the galaxy angular size.  If not, we can ignore the correction.
- Calculate a new grid roughly centered on the localization but shifted from the galaxy grid.  That is, the localization grid is like an extension of the galaxy grid.  
- I have put a hook in the code to allow for a correction.  If b < phi_cand, you should calculate the correction and apply it to the grid.
- Log your work in the "Logs" section below.

4. While your _local_lwx_correction helper works well, it is not exactly what I had in mind.  Let us develop an alternative approach which won't be as accurate but should be much faster.  Here is what the helper function should do:

- Calculate L_wx approximately centered on the localization.  
- The grid should be large enough to cover most of the ellipse (4 sigma along the major axis) 
- The grid is spaced exactly the same as the galaxy grid (step_size_phi)
- The grid should be aligned to the galaxy grid, i.e. only shifted by an integer number of cells in x and y
- For phi_cand much greater than b, the total L_wx will be significantly less than 1 and we will correct our calculation by this factor
- We will name this helper _Lwx_correction() 
- Be careful not to generate any arrays larger than approximately 5000x5000 cells.

The main point is that we will under-resolve the localization in our calculation.  But we can calculate the correction factor and apply it to the grid.  I have found that this is accurate to ~1% for the cases I have tested.

Consider the above and Log a plan in the "Logs" section below.  Do not generate any new code yet.

5. What you describe in the logs looks correct.  Please:

- proeed to generate the code.  
- Do use a square window for the correction.
- Log your work in the "Logs" section below.


## Profiling

1. Please make the following modification to profiling.py:

- Add an additional test for the local method which uses an ellipse localization with a relatively large axis ratio.  Large enough to require an approximately 1000x1000 grid for the L_wx correction.
- Include this curve in the figure for the local method.

2. The purple curve in the local figure says circular but no correction.  Please:

- check whether there is a correction.  
- update the figure as needed
- Log your work in the "Logs" section below.

3. Add one more test to profiling.py for the local method:

- Use a very small localization in both dimenstions (a=0.1", b=0.1")
- Include galaxy sizes from 0.1" to 20"
- If necessary, have the test_local.py module include an accuracy test for a small localization too
- Log your work in the "Logs" section below.

## Docs


1. Examine the docs in the astropath/docs directory.  We need to update the documentation to reflect the new code.  Please:

- Update the documentation to reflect the new code for the local method
- Note that it does not require numba (for now)
- Comment on the profiling.py module too
- Log your work in the "Logs" section below.

2. Please add the results on accuracy and profiling to the docs/performance.rst file.  And:

- Include the 2 profiling figures (re-generate them first)
- Tables summarizing the accuracy and profiling results

## Prompts

1. Read this doc.  Proceed with the 1st item under Testing.
2. Read this doc.  Proceed with the 1st item under Development.
3. Read this doc.  Proceed with the 2nd item under Development.
4. Read this doc.  Proceed with the 3rd item under Development.
5. Read this doc.  Proceed with the 4th item under Development.
6. Read this doc.  Proceed with the 5th item under Development.
7. Read this doc.  Proceed with the 1st item under Docs
8. Read this doc.  Proceed with the 1st item under Profiling
9. Read this doc.  Proceed with the 2nd item under Profiling
10. Read this doc.  Proceed with the 3rd item under Profiling
11. Read this doc.  Proceed with the 2nd item under Docs

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

### 2026-06-11 (Small-localization correction, centered on the loc)

**Done** (Development item 3).  The default step is now 0.05 (relative
to phi).  That is fine when the localization is comparable to / larger
than the galaxy, but when the minor axis `b < phi` the galaxy-centered
grid (spacing `phi*step`) under-resolves the sharp localization and the
raw sum carries the full O(step) bias (~-0.8% at step 0.05, ~-1.7% at
0.1).  I filled the existing `if b < phi_cand:` hook with a correction.

**What the correction does** (new helper
`bayesian._local_lwx_correction`):
- Builds a SEPARATE grid CENTERED ON THE LOCALIZATION (the transient),
  not the galaxy -- the galaxy grid shifted/extended onto the
  localization, in the same flat-sky tangent plane.
- The grid is tied to the ELLIPSE, not to `step_size`: half-width
  `6a` (covers the major axis so `sum(L)` captures the full,
  unit-normalized L) and spacing `b/4` (resolves the minor axis).  It is
  therefore INDEPENDENT of `step_size`, so coarsening the galaxy step
  no longer degrades these cases.  `ngrid` is capped at 4000 for extreme
  axis ratios.
- Returns `sum(L*p)/sum(L)` -- the L-weighted average of p(w|O_i).  The
  grid spacing cancels in the ratio and L integrates to 1, so this is
  the localization-normalized p(x|O_i), exactly the quantity
  `px_Oi_fixedgrid` produces with `correction='L_wx'`.  This is "similar
  to the fixed method", but centered on the localization as required.
- Only the `eellipse` path uses it (per item 2's scope).  For `b >=
  phi` nothing changes; the generic healpix/wcs path is untouched.

**Why a separate, loc-centered grid (the subtlety):** the ratio's
denominator `sum(L)` must capture the FULL localization (so it equals
its analytic norm of 1); a galaxy-centered grid of half-width `phi*max`
fails this when `a > phi*max` (e.g. the a=10",b=0.2" ellipse), giving
`sum(L)` over only a sliver and a ratio ~4x too large.  Centering on the
localization with half-width `6a` fixes it.

**Empirical design** (vs the fine fixed-grid reference, at step 0.05):
| case (b<phi)            | raw reldiff | corrected | ngrid |
|-------------------------|-------------|-----------|-------|
| large gal10 / loc1      | -8.5e-3     | -2.7e-5   | 48    |
| gal10 / loc0.3          | -1.2e-2     | -1.3e-7*  | 48    |
| ellipse a10",b0.2"      | -8.4e-3     | -8.0e-5   | 2400  |
(*with the chosen b/4 spacing; a step-tied spacing under-resolved this
case, which is why the correction grid is tied to the ellipse.)  A
galaxy-grid L_wx ratio was also tried but blows up (+3.2) on the ellipse
because the galaxy box doesn't cover the major axis -- hence the
loc-centered grid.

**Verification:** `pytest astropath/tests/tests_local.py -s` -> 11
passed.  Added `test_local_correction_coarse_step` (parametrized over
the two `b<phi` cases x step_size in {0.05, 0.10}): asserts agreement to
`rtol=2e-3`, far tighter than the generic 1% -- and the corrected value
is identical at 0.05 and 0.10, confirming step-independence.
`test_bayesian` and `test_path` still pass.

**Unrelated note:** `astropath/run.py:231` has a stray
`embed(header='run.py:231')` (in the committed `wip` state, hit
unconditionally for both local and fixed modes) that hangs
`test_run.py::test_run_on_dict_*` under pytest with
"OSError: reading from stdin while output is captured".  I left it in
place (it is your WIP debugging, not part of this task); remove that
line to make those two tests pass.

**Note / possible follow-up:** for very high axis-ratio ellipses
(`a >> b`) the correction grid is `~6a x 6a` at spacing `b/4`, i.e.
`~24 a/b` pixels per side (2400 for the test ellipse; the 4000 cap can
bite for, say, a=60",b=0.2").  An anisotropic grid aligned to the
ellipse axes (fine along the minor, coarse along the major) would cut
this by ~100x; I kept the uniform grid for fidelity to "an extension of
the galaxy grid" and simplicity.  Flag if you want the anisotropic
version.

### 2026-06-11 (PLAN: faster aligned _Lwx_correction; no code yet)

**Task** (Development item 4): replace the accurate-but-slow
`_local_lwx_correction` with a faster, ~1%-accurate `_Lwx_correction()`
that DELIBERATELY under-resolves the localization on the (coarse) galaxy
grid and divides the raw result by a correction FACTOR.  The hook is
already stubbed by you:
```
if b < phi_cand:
    L_wx_correction = _Lwx_correction(arguments)
else:
    L_wx_correction = 1.0
...
grid_p = L_wx * p_wOi / L_wx_correction          # scalar divide
p_xOis.append(np.sum(grid_p) * step_size_phi**2)
```
So `_Lwx_correction()` must return a SCALAR factor; `p(x|O_i)` becomes
`[sum(L*p)*dA] / factor`.

**What the factor is.**  `factor = sum_over(loc grid) L_wx * dA`, i.e.
the DISCRETE integral of L ("total L_wx") evaluated on a coarse grid
centered on the localization.  When `phi >> b` the galaxy-grid spacing
`step_size_phi` under-samples the sharp ellipse, so this discrete
integral comes out `< 1` (your observation).  Dividing the raw
`sum(L*p)*dA` by it renormalizes L back to 1.

**Why this works to ~1% even though L is under-resolved -- the key
idea is GRID ALIGNMENT.**  The raw numerator `sum_gal L*p*dA` samples L
on the galaxy lattice at some fixed sub-cell phase relative to the
ellipse.  If the correction grid is the SAME lattice (same spacing,
shifted by an integer number of cells), it samples L at the SAME
sub-cell phase.  Under-resolution then multiplies BOTH the numerator's
L-weighting and the denominator `sum_loc L*dA` by the same (roughly
uniform) aliasing factor, so it cancels in the ratio:

    [sum_gal L*p*dA] / [sum_loc L*dA]
        ~= integral(L*p) / integral(L)  =  integral(L*p)  =  p(x|O_i)

(the last step uses integral(L)=1).  The residual ~1% comes from the
aliasing factor not being perfectly uniform, the 4-sigma truncation,
and p varying across the product region.  Alignment is ESSENTIAL: an
unaligned (differently-phased) correction grid would not cancel and the
error would be large and noisy.  This is the crucial difference from
`_local_lwx_correction`, which instead RESOLVED L on a fine (b/4) grid
-- accurate (~1e-5) but slow.

**Planned `_Lwx_correction()` -- spec.**
- Purpose: return the scalar `sum(L_wx)*dA` over a coarse,
  galaxy-grid-aligned grid centered on the localization, so the caller
  can divide the (under-resolved) raw p(x|O_i) by it.
- Inputs (all plain floats/arrays, numba-friendly):
  * `E0, N0` -- flat-sky offsets of the galaxy from the localization
    center (arcsec); already computed in the loop.
  * `a, b, cos_dth, sin_dth` -- ellipse semi-axes and the shared
    rotation into the ellipse frame.
  * the galaxy-grid lattice description so the correction grid can be
    built ON the same lattice: pass `box_hwidth` and the shared
    normalized axis `u` (=`np.linspace(-1,1,ngrid)`), from which
    `h = box_hwidth*(u[1]-u[0])` is the EXACT galaxy-grid spacing
    (~`step_size_phi`).  (Using the exact `h`, not the nominal
    `step_size_phi`, keeps the two grids truly aligned.)
  * `step_size_phi` -- only for `dA = step_size_phi**2` (kept identical
    to the caller's `dA` so the scalar divide is consistent; any common
    dA cancels in the ratio anyway).
- Grid construction (the three constraints you listed):
  1. SAME spacing as the galaxy grid: spacing `= h`.
  2. ALIGNED / integer-cell shift: snap the localization center onto the
     galaxy lattice -> `iE = round(E0/h)`, `iN = round(N0/h)`; the 1D
     correction coordinates (offsets from the galaxy) are
     `(iE + arange(-m, m+1)) * h` and likewise for N, where
     `m = ceil(4*a / h)`.  These points lie exactly on the galaxy
     lattice, just translated by `(iE, iN)` cells.
  3. COVER 4 sigma of the major axis: half-width `4*a` -> the
     `m = ceil(4a/h)` above.  (A square of half-width 4a covers 4-sigma
     in every direction, so PA need not be considered for sizing.)
- Evaluate L on this grid: offset of each correction-grid point from the
  LOCALIZATION center is `(cE - E0, cN - N0)`; rotate into the ellipse
  frame (`x_box = E*cos_dth + N*sin_dth`, `y_box = N*cos_dth -
  E*sin_dth`) and form the same 2D Gaussian
  `exp(-x_box^2/2a^2)*exp(-y_box^2/2b^2)/(2*pi*a*b)` already used in the
  main loop.  Reuse, don't duplicate, that expression where practical.
- Return `np.sum(L_corr) * step_size_phi**2`  (the scalar "total L_wx").
- Guard: if `b >= phi` the caller never calls this (factor=1.0).  Cap
  `m` (e.g. `2*m+1 <= ~4000`) so extreme axis ratios stay bounded.

**Cost vs `_local_lwx_correction`.**  New spacing is `h ~ phi*step`
(coarse) instead of `b/4` (fine).  Correction-grid side length
`~ 8a/h = 8a/(phi*step)` vs the old `~ 48a/b`.  Ratio
`new/old ~ b/(6*phi*step)`, i.e. for `phi >> b` (the regime that
triggers the correction) the new grid is many-fold smaller -- e.g.
phi=10",b=0.2",step=0.05 -> ~15x fewer points; and it sums only L (no
p, no product) per point.  For low-phi high-aspect ellipses
(phi~b, e.g. the a=10",b=0.2",phi=0.5" test case) the two are
comparable, since there `h` is already small.

**Integration / cleanup.**
- Implement `_Lwx_correction()` with the spec above; wire its real
  arguments into your stub (`E0, N0, a, b, cos_dth, sin_dth,
  box_hwidth, u, step_size_phi`).
- Delete `_local_lwx_correction` (superseded) -- per your earlier note
  we need not preserve old algorithms; `tests_local.py` is the anchor.

**Test impact (to handle in the coding step).**
- `test_local_correction_coarse_step` currently asserts `rtol=2e-3`
  (true of the fine helper).  The new approach is ~1%, so relax that to
  ~`1.5e-2` (and drop the "identical at 0.05 vs 0.10" expectation -- the
  factor now depends on the grid, hence weakly on step).  I will measure
  the actual reldiff for the two `b<phi` cases at steps {0.05,0.10} and
  set the tolerance just above the worst observed.
- `test_local_matches_fine_fixedgrid` runs at step 0.02, where the
  galaxy grid already resolves L for these cases, so its `rtol=1e-2`
  should still pass; I will confirm.

**numba note.**  `_Lwx_correction` is a small Gaussian sum over a
coarse aligned grid -- trivially `@njit`-able later, and (like
`px_Oi_numba`) it takes only scalars/arrays, no astropy.

**Open question (will proceed with the stated default unless you say
otherwise):** "4 sigma along the major axis" -> I read this as
half-width `= 4*a` (with `a` the major-axis sigma).  If you instead
meant 4 sigma in each axis separately (an anisotropic/elliptical
window), say so; the square 4a window is simpler, stays aligned, and
covers the ellipse, so I will use it.

### 2026-06-11 (Implemented fast aligned _Lwx_correction)

**Done** (Development item 5).  Replaced `_local_lwx_correction` with
the faster `bayesian._Lwx_correction()` and wired it into your stubbed
hook.  Square 4-sigma window, as you confirmed.

**`_Lwx_correction(E0, N0, a, b, cos_dth, sin_dth, box_hwidth, ngrid,
step_size_phi, max_side=5000)`** returns a SCALAR "total L_wx" =
`sum(L_wx)*dA`.  The caller does `grid_p = L_wx*p_wOi /
L_wx_correction`, i.e. `p(x|O_i) = raw / factor`.  The grid:
- is centered on the localization (transient), square half-width `4*a`
  (>= 4 sigma of the major axis);
- has EXACTLY the galaxy-grid spacing `h = 2*box_hwidth/(ngrid-1)` and
  lies ON the galaxy lattice (snapped by an integer cell shift `kE, kN`
  via `round((box_hwidth - E0)/h)`), so L is sampled at the same
  sub-cell phase as the raw sum -> the under-resolution aliasing cancels
  in the `raw/factor` ratio;
- uses `dA = step_size_phi**2` to match the caller's raw sum exactly
  (so both the aliasing AND the step_size_phi-vs-h convention cancel).

**Memory guard (your <=5000x5000 constraint).**  The window has
`~8*a/h` cells per side; this only exceeds `max_side=5000` when `h` is
very small, which is exactly the regime where the galaxy grid ALREADY
resolves L and the raw value needs no correction.  So when
`2m+1 > max_side` the helper returns `1.0` (skip) -- provably safe and
never allocates a large array.  (Found empirically while validating:
the ellipse at step 0.02 wants a 7987-cell window -> skipped -> raw is
already -0.33%.)  No striding/coarsening: that would break the
aliasing/dA cancellation, so skip-when-resolved is the right lever.

**Accuracy measured vs the fine fixed grid** (reldiff; raw = no
correction):

| case (b<phi)          | step | raw     | corrected | window |
|-----------------------|------|---------|-----------|--------|
| gal10", loc 1"        | 0.05 | -0.85%  | -1.7e-4   | 17     |
| gal10", loc 1"        | 0.10 | -1.85%  | -1.9e-3   | 9      |
| gal10", loc 1"        | 0.02 | -0.33%  | +7.0e-5   | 41     |
| ellipse a10",b0.2"    | 0.05 | -0.84%  | -4.1e-6   | 3189   |
| ellipse a10",b0.2"    | 0.10 | -1.66%  | +3.4e-5   | 1589   |
| ellipse a10",b0.2"    | 0.02 | -0.33%  | -0.33%(skip)| --    |
| severe UR phi2 b0.05  | 0.05 | -2.49%  | -1.8e-3   | 9      |
| severe UR phi2 b0.05  | 0.10 | -60.1%  | -2.1%     | 5      |

So it removes the bulk of the under-resolution bias (dramatically so in
the severe cases), staying ~0.1-0.2% at the default step 0.05 -- better
than the ~1% you found.  (The lone outlier, severe under-resolution at
the very coarse step 0.10 with only 5 window cells, lands ~2%.)

**Cost.**  Common large-galaxy/small-loc case: window ~17-41 cells,
negligible.  px_Oi_local at step 0.05 with the correction firing on
ALL 100 candidates (a=1", phi=2-8"): 0.74 ms/candidate.  Much cheaper
than the old fine-grid helper (which used a b/4 spacing -> e.g. 48-cell
side even for gal10"/loc1" and 2400 for the ellipse).

**Tests.**  `tests_local.py` -> 11 passed.
- `test_local_matches_fine_fixedgrid` (step 0.02, rtol 1e-2): the two
  b<phi cases pass -- gal10"/loc1" via the correction (+7e-5), the
  ellipse via the skip path (-3.3e-3, raw).  No large array built.
- `test_local_correction_coarse_step`: relaxed rtol 2e-3 -> 5e-3 (the
  new method is ~1%-class, not the old ~1e-5) and dropped the
  step-independence wording (the factor now depends on the galaxy grid).
  Worst observed: gal10"/loc1" at step 0.10 = -1.9e-3.
- Removed the obsolete `_local_lwx_correction`.
- `test_bayesian`, `test_path` still pass.

**numba note.**  `_Lwx_correction` is a small Gaussian sum on a coarse
aligned grid (scalars + arrays, no astropy) -- ready for an `@njit`
kernel in the numba step.

### 2026-06-11 (Docs: documented the local method in performance.rst)

**Done** (Docs item 1).  The docs live in the top-level `docs/`
directory (not `astropath/docs`).  `docs/performance.rst` already
documented the numpy-accelerated core, the optional numba path, and the
`profiling` module (from the earlier fixed-grid work), so it was the
natural home for the local-method update.

**Edits to `docs/performance.rst`:**
- Added a third bullet to *numpy-accelerated core* noting that
  `px_Oi_local` was rewritten in pure numpy for the `eellipse` type
  (coords pre-extracted once, `L(w-x)` evaluated flat-sky in the tangent
  plane, no astropy in the per-candidate loop; healpix/wcs still use
  `calc_LWx`), with a cross-reference to the new section.
- New *Local-grid posterior* section: what `px_Oi_local` is for (large
  localizations / healpix), the once-built normalized grid
  (`ngrid = 2*max/step_size`), the flat-sky `L(w-x)` matching `calc_LWx`
  to ~1e-5, and that `step_size` is *relative* to the galaxy (default
  0.05).  Includes a prominent **note that it is pure numpy and does
  NOT require/use numba** ("for now"); numba remains fixed-grid only.
- New *Small-localization correction* subsection describing the
  `b < phi` case and `_Lwx_correction`: the localization-centered,
  galaxy-aligned "total L_wx" factor, why the aliasing cancels (~1%),
  the analogy to `px_Oi_fixedgrid`'s `correction='L_wx'`, and that the
  correction grid is bounded (skipped above ~5000 cells/side, which only
  happens when no correction is needed) so it never allocates a large
  array.
- Updated the *Profiling module* section to say the sweep now profiles
  BOTH methods: the existing fixed-grid (+numba) curve
  (`profiling_timing.png`) and the new `run_profiling_local` /
  `plot_local_results` local-grid curve (`profiling_local_timing.png`),
  with both tables printed.  Added the import-and-call example for the
  local profiling.

**Verification:** installed `sphinx`/`sphinx_rtd_theme`/`nbsphinx` in
astro14 and ran `python -m sphinx -b html . _build/html`.  Build
succeeded; the 11 warnings are all pre-existing and unrelated (missing
`nb/*` notebooks, `_static`, `language=None`, the `run.py` docstring
list, `chime.rst` not in a toctree).  No warning references
`performance.rst`.  Confirmed in the generated `performance.html` that
the new *Local-grid posterior* / *Small-localization correction*
sections, the no-numba note, the profiling mentions, and the
`Local-grid posterior`_ cross-reference (-> `#local-grid-posterior`) all
render.

### 2026-06-11 (Profiling: ellipse scenario exercising _Lwx_correction)

**Done** (Profiling item 1).  Added a high-axis-ratio ellipse scenario
to the local-method profiling so the `_Lwx_correction` path is timed on
a ~1000x1000 grid, and plotted its curve alongside the existing one.

**Edits to `astropath/profiling.py`:**
- New `ellipse_setup()`: a long thin error ellipse (a=12.5", b=0.2",
  axis ratio ~62) with 50 candidate galaxies of size 1.5"-2.5" (all
  `> b`, so the correction fires for every candidate).  The galaxy
  sizes are tuned so the correction grid is ~1000 cells per side at the
  default step 0.05, and stays under the ~5000-cell skip threshold
  across the whole step sweep (so it never short-circuits).
- `run_profiling_local()` now times BOTH scenarios against the same
  per-candidate galaxy grid (the figure x-axis): the existing circular
  case (`b >= phi`, no correction) and the new ellipse case
  (`b < phi`, correction fires).  Added columns
  `px_Oi_local_ellipse_s` and `corr_ngrid` (a representative
  correction-grid side, ~`8a/h`, printed per step), so the table shows
  the correction grid growing 97 -> 1777 across the sweep
  (997 at the default 0.05).
- `plot_local_results()` plots the ellipse curve (orange squares,
  "ellipse + L_wx corr.") on top of the circular curve (purple
  diamonds), keeping the 10 s reference line.
- `main()`'s local table now reports `corr_ngrid`, `circular_ms`, and
  `ellipse_ms`.

**Measured** (50 candidates; circular vs ellipse+correction):

| step  | galaxy grid | corr grid | circular | ellipse(+corr) |
|-------|-------------|-----------|----------|----------------|
| 0.50  | 24x24       | 97        | 1.0 ms   | 6.7 ms         |
| 0.25  | 48x48       | 197       | 1.9 ms   | 32 ms          |
| 0.10  | 120x120     | 497       | 7.9 ms   | 361 ms         |
| 0.05  | 240x240     | 997       | 32 ms    | 1763 ms        |
| 0.025 | 480x480     | 1777      | 157 ms   | 7059 ms        |

So the L_wx correction dominates the cost for high-axis-ratio ellipses
(the correction grid `~8a/h` grows faster than the galaxy grid
`2*max/step`): at step 0.05 the per-candidate correction is a ~1000x1000
Gaussian sum, ~35 ms/candidate vs ~0.6 ms for the un-corrected circular
case.  This is the intended demonstration that the cheap-but-not-free
correction is the bottleneck precisely for long thin localizations --
useful context for a future numba kernel.

**Verification:** ran `run_profiling_local` and `plot_local_results`;
confirmed `corr_ngrid=997` at step 0.05, that both data curves render
(labels "px_Oi_local (circular, no corr.)" and "px_Oi_local (ellipse +
L_wx corr.)"), and that the largest correction array stays well under
5000x5000 (max ~2663 per side at step 0.025).

### 2026-06-11 (Profiling: verified the circular curve has no correction)

**Done** (Profiling item 2).  Checked whether the purple (circular)
curve in the local figure actually applies the `_Lwx_correction`.

**Verification:** wrapped `bayesian._Lwx_correction` with a call
counter and ran both profiling scenarios at step 0.05:
- circular (`default_setup`, a=b=5"): the correction fired **0 times**
  -- `b < phi` is False for every candidate (b=5" vs phi 0.2"-2.5"), so
  the `else: L_wx_correction = 1.0` branch is always taken.  The
  "no correction" label is therefore correct.
- ellipse (`ellipse_setup`, b=0.2"): the correction fired **50 times**
  (every candidate), with factor ~0.9916.

So the purple curve genuinely has no correction; the difference between
the two curves IS the correction overhead.

**Figure update:** sharpened both legend labels to state the governing
condition explicitly, so the figure is self-documenting:
- ``px_Oi_local (circular, b>=phi: no correction)``
- ``px_Oi_local (ellipse, b<phi: L_wx correction)``
and added comments in `plot_local_results` noting the circular case has
`b >= phi` (0 invocations, verified) while the ellipse case fires the
correction on a ~1000x1000 grid.

### 2026-06-11 (Profiling: tiny-localization scenario + NaN fix)

**Done** (Profiling item 3).  Added a third local-method profiling
scenario and an accuracy test; the scenario also surfaced (and I fixed)
a 0/0 NaN in `px_Oi_local`.

**`ellipse_setup` companion -> new `small_loc_setup()`:** a very small
circular localization (a=b=0.1") with galaxy sizes 0.1"-20".  For every
galaxy `> b` the `_Lwx_correction` fires, but because the major axis is
tiny the correction grid stays small (the opposite of `ellipse_setup`'s
~1000-cell grid).  This stresses the deeply under-resolved regime
(galaxy spacing `phi*step` can be many times `b`).

**`run_profiling_local` now times THREE scenarios** against the same
per-candidate galaxy grid (figure x-axis): circular (no correction),
ellipse (large correction grid), and small-loc (tiny correction grid).
Added column `px_Oi_local_smallloc_s`; `plot_local_results` plots the
new curve (green triangles) and `main()`'s table gains `smallloc_ms`.
Full sweep (50 candidates):

| step  | galaxy grid | ellipse corr | circular | ellipse | small-loc |
|-------|-------------|--------------|----------|---------|-----------|
| 0.50  | 24x24       | 97           | 1.0 ms   | 6.8 ms  | 2.3 ms    |
| 0.25  | 48x48       | 197          | 1.8 ms   | 32 ms   | 3.3 ms    |
| 0.10  | 120x120     | 497          | 7.9 ms   | 354 ms  | 11 ms     |
| 0.05  | 240x240     | 997          | 32 ms    | 1822 ms | 40 ms     |
| 0.025 | 480x480     | 1997         | 156 ms   | 9721 ms | 184 ms    |

So the correction cost scales with the localization MAJOR axis (a): a
long thin ellipse (a=12.5") is ~50x the no-correction cost, while a tiny
round localization (a=0.1") adds only a small overhead.

**Bug found + fixed (0/0 -> NaN).**  The small-loc scenario at coarse
steps produced NaN: when the grid spacing `phi*step` is many times
`b` (e.g. step 0.5, phi up to 20", spacing up to 10" vs b=0.1"), the
0.1" localization falls entirely between grid cells -- the aligned
`_Lwx_correction` grid catches no flux (factor underflows to 0) and the
raw sum is 0 too, so `raw/factor` was `0/0`.  This is a real pipeline
risk (tiny localization + coarse step), so I guarded `px_Oi_local`: when
the correction factor is not `> 0`, it falls back to the exact
delta-function limit `p(x|O_i) = pw_Oi(theta_offset)` (correct as
`b -> 0`, since L integrates to 1).  After the fix the small-loc
scenario is finite at every step (verified 0.5 -> 0.025, no warnings).

**Accuracy test added** to `tests_local.py`
(`test_local_small_localization`, 3 cases): a=b=0.1", phi in
{0.3, 0.6, 1.0}", offset 0.3", compared to the fine fixed grid at the
default step 0.05.  phi is kept <= 1" so the reference stays well under
5000 cells per side.  Result: reldiff ~2e-5 to 7e-5 -- the correction
is essentially exact for very small localizations.  (The earlier ~-2%
seen against a delta-limit reference was that reference's own bias, not
the method.)  Full suite: `tests_local.py` 14 passed; with
`test_bayesian`/`test_path`, 19 passed.

**Note:** the largest array anywhere stays under 5000x5000 -- the
small-loc correction grids are tiny (small a), the ellipse tops out at
~2000 (step 0.025), and the NaN guard removed the only pathological
divide.
