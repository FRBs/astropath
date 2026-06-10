# Speeding up the PATH calculations

## Goals

This document is a collection of ideas and prompts for speeding up the PATH calculations. 

## Coding

Here are guidelines for coding: 

- Use Python where possible
- Add inline comments to explain the effort
- Reuse existing code when possible
- Place import statements at the top of the file.
- Include a description of inputs/outputs in the doc string of all methods
- Use lines of code that are less than 80 characters wide

## Running

If you need to run python, use the "astro14" environment.  

## Development

### calc_LWx

We wish to speed up the calculation of L(w-x) for a given localization. 

1. Speed up the code in localization.py:calc_LWx that calculates the grid of L(w-x) values for a given localization with the eellipse type.  Please:

- Examine the code
- Identify the bottlenecks
- Figure out how to replace all astropy calculations 
- Log your work in the "Logs" section below.
- Make recommendations 
- Do not generate any new code yet
- If you have any questions, ask them in the Logs and I will clarify at my next prompt

2. Write new code to replace the astropy calculations with numpy calculations. Please:

- Replace the astropy calculations with numpy calculations as described in the Logs.
- Keep the existing astropy code as a test case.
- Stick with calc_Lwx for now;  we will revisit px_Oi_local later.
- Generate a new test case to verify the new code.  Make a test_speed_up.py module and model it after the other modules in astropath/tests.
- Have the tests also compare and report the speed up
- Include the timing for a grid that is 7200x7200 pixels
- Log your work in the "Logs" section below.

### px_Oi_fixedgrid

We wish to speed up the calculation of p(x|O_i) for a given localization using the fixed grid method.  

1. Similar to the calc_LWx bottlenecks, we need to identify the bottlenecks and replace the astropy calculations with numpy calculations.  Please:

- Examine the code
- Identify the bottlenecks
- Figure out how to replace all astropy calculations 
- If it would help to store ra, dec in radians and not deg consider it
- Log your work in the "Logs" section below.
- Make recommendations 
- Do not generate any new code yet
- If you have any questions, ask them in the Logs and I will clarify at my next prompt

2. Write new code to replace the astropy calculations with numpy calculations. Please:

- Replace the astropy calculations with numpy calculations as described in the Logs.
- Generate new tests and add them to the test_speed_up.py module.
- Only implement items 1 and 2 in the Logs; no vectorization or other optimizations yet
- We may try to speed things up further with numba.  Consider that in your programming
- Log your work in the "Logs" section below.

## Profiling

1. Generate a module to profile the code.  Call it profiling.py and place it in the astropath directory.  At this stage, we wish to profile the code for the calc_LWx and px_Oi_fixedgrid methods.  Please:

- Model the code after the Profiling.ipynb notebook in the astropath/calculations/step_size directory
- Run on a range of grid sizes by varying the step_size
- Generate a table of the timing results.  Print these to the screen
- Generate a figure of the timing results 
- Log your work in the "Logs" section below.

2. Modify the profiling.py tests to:

- Run on N=50 candidates instead of only 2.  Have the galaxies have a range of anuglar sizes and locations
- 

## Prompts

1. Read this doc.  Proceed with the 1st item under Development/calc_LWx
2. Read this doc.  Proceed with the 2nd item under Development/calc_LWx
3. Read this doc.  Proceed with the 1st item under Development/px_Oi_fixedgrid
4. Read this doc.  Proceed with the 2nd item under Development/px_Oi_fixedgrid
5. Read this doc.  Proceed with the 1st item under Development/Profiling


## Logging

The "Logs" section will record Claude's work.  Please use the following format:

### <Date> (Short summary of the work)

<Detailed description of the work and what you learned>

### <Date> (Short summary of the work)

<Detailed description of the work and what you learned>

...


## Logs

### 2026-06-09 (Analysis of calc_LWx eellipse bottlenecks)

**Code examined:** `astropath/localization.py:calc_LWx`, the `eellipse`
branch (lines 63-81), plus its two callers in
`astropath/bayesian.py`: `px_Oi_fixedgrid` (line 122, called once per
analysis) and `px_Oi_local` (line 209, called **once per candidate**
inside the galaxy loop).

**What the eellipse branch does:**
1. `SkyCoord(ra=ra, dec=dec, unit='deg')` — builds a SkyCoord from the
   whole grid (lines 69-70).
2. `center.separation(coord)` — Vincenty angular separation (line 72).
3. `center.position_angle(coord)` — position angle E of N (line 73).
4. Rotates into the ellipse frame and evaluates a 2D Gaussian (74-81).

**Bottleneck identified:** steps 1-3 are the cost. They are pure
astropy and dominate the branch. The arithmetic in steps 1-3 is
trivial; the expense is astropy's per-call machinery — SkyCoord object
construction, `Quantity`/unit wrapping, and `.to('arcsec')` /
`.to('deg')` conversions. I confirmed this with two benchmarks
(grid mimicking `box_hwidth=10"`, `step_size=0.1"`):

| grid        | pts   | astropy   | numpy     | speedup |
|-------------|-------|-----------|-----------|---------|
| 20x20       | 400   | 0.342 ms  | 0.018 ms  | ~19x    |
| 60x60       | 3600  | 0.509 ms  | 0.112 ms  | ~4.5x   |
| 200x200     | 40000 | 3.13 ms   | 1.13 ms   | ~2.8x   |

The speedup grows as the grid shrinks because astropy's overhead is
largely fixed per call. This matters most in `px_Oi_local`, which
calls `calc_LWx` once per candidate on small per-galaxy grids — the
regime where the win is ~19x.

**How to replace all astropy calculations:** astropy's `separation` and
`position_angle` are just closed-form spherical formulas
(`astropy.coordinates.angle_utilities`). They can be reproduced exactly
in numpy. Working from the center (ra0, dec0) and grid (ra, dec), all
in radians:

```
dlon = ra - ra0
sdlon, cdlon = sin(dlon), cos(dlon)
sl1, cl1 = sin(dec0), cos(dec0)
sl2, cl2 = sin(dec), cos(dec)
# separation (Vincenty), radians:
sep = arctan2(hypot(cl2*sdlon, cl1*sl2 - sl1*cl2*cdlon),
              sl1*sl2 + cl1*cl2*cdlon)
# position angle E of N, radians:
pa  = arctan2(sdlon*cl2, sl2*cl1 - cl2*sl1*cdlon)
```

I verified these against astropy on the 200x200 grid: max separation
diff ~7e-11 arcsec, max PA diff ~1e-8 deg — i.e. bit-for-bit
equivalent for our purposes. `sep` then goes to arcsec via
`degrees(sep)*3600`, `pa` to degrees via `degrees(pa)` (no need to
wrap to [0,360) since the result is fed straight into `sin`/`cos`).

The center RA/Dec come from `localiz['center_coord'].ra.deg` /
`.dec.deg` (read once, cheap). The `coord.equinox =
localiz['center_coord'].equinox` line (70) is a no-op for the
math — ICRS separation/PA don't depend on equinox — and can be
dropped in the numpy path.

**Recommendations (no code written yet):**
1. Replace the SkyCoord/separation/position_angle block in the
   `eellipse` branch with the numpy formulas above. Identical results,
   ~3-19x faster depending on grid size.
2. Bigger structural win: in `px_Oi_local` the localization mapping is
   rebuilt per candidate. `px_Oi_fixedgrid` already computes `L_wx`
   once (good). Worth revisiting `px_Oi_local` after step 1.
3. Keep the spherical (not flat-sky) formula so numerical results are
   unchanged; only the implementation changes.

**Question for next prompt:** Should I prioritize (a) just the in-place
numpy swap inside `calc_LWx`'s eellipse branch (lowest risk, drop-in,
identical output), or also (b) vectorize/cache across the candidate
loop in `px_Oi_local`? And is exact equivalence with astropy required,
or is the existing flat-sky approximation already accepted such that I
could simplify further?

### 2026-06-09 (Implemented numpy calc_LWx + speed-up tests)

**Code changed:** `astropath/localization.py:calc_LWx`, eellipse branch.
Replaced the astropy block (`SkyCoord` build, `.separation()`,
`.position_angle()`, and `units`-based rotation) with the pure-numpy
Vincenty separation + position-angle formulas worked out in the prior
log entry. The rest of the branch (rotation into the ellipse frame, the
2D-Gaussian evaluation) is unchanged, so output is numerically
identical. The `coord.equinox = ...` line was dropped (no-op for ICRS
separation/PA). Note: the `from astropy import units` import in
localization.py is now unused (a Hint-level lint); left in place since
it is harmless and `units` is a common re-export — can remove if
desired.

**Tests added:** `astropath/tests/test_speed_up.py`, modeled on
`tests/test_localization.py`. It freezes the original astropy eellipse
code as `_calc_LWx_eellipse_astropy()` (the requested "keep the astropy
code as a test case") and:
- `test_eellipse_matches_astropy` — 200x200 grid, asserts
  `np.allclose(rtol=1e-8, atol=1e-12)` against the astropy reference.
- `test_eellipse_matches_astropy_offcenter` — southern, high-PA case
  (dec=-50.8, theta=110) to exercise the geometry.
- `test_eellipse_speed_up[200|1000|7200]` — times both implementations,
  re-checks correctness, prints the speed-up, and asserts numpy is not
  slower. Includes the requested **7200x7200** grid.

**Measured results** (astro14 env, best-of-3, run with `pytest -s`):

| grid        | pts        | astropy    | numpy      | speed-up |
|-------------|------------|------------|------------|----------|
| 200x200     | 40,000     | 2.95 ms    | 2.39 ms    | 1.2x     |
| 1000x1000   | 1,000,000  | 141.3 ms   | 85.0 ms    | 1.7x     |
| 7200x7200   | 51,840,000 | 11,984.9 ms| 6,092.0 ms | 2.0x     |

For very large grids the win plateaus at ~2x (memory-bandwidth-bound on
the shared trig terms, no longer astropy-overhead-bound). The big-grid
gains are smaller than the small-grid microbenchmark in the previous
entry (~19x at 20x20) precisely because the per-call astropy overhead
is fixed — it dominates only when the grids are small, i.e. the
`px_Oi_local` per-candidate regime (to be revisited later).

**Verification:** all 5 new tests pass; the existing
`test_path.py`, `test_localization.py`, and `test_bayesian.py` suites
(7 tests) still pass, confirming the eellipse path is unchanged
downstream.

### 2026-06-09 (Analysis of px_Oi_fixedgrid bottlenecks)

**Code examined:** `astropath/bayesian.py:px_Oi_fixedgrid`
(lines 67-153), plus the `pw_Oi` it calls (15-64).

**Astropy usage in the function:**
1. Line 105 — `localiz['center_coord'].equinox = cand_coords[0].equinox`
   (one-time; sets an attribute).
2. Lines 119-120 — `localiz['center_coord'].ra.deg`,
   `np.cos(localiz['center_coord'].dec).value` to build the grid
   (one-time, scalar; cheap).
3. Line 122 — `calc_LWx` (now pure numpy for eellipse; one call).
4. **Lines 126-130 (the per-candidate loop)** —
   `for cand_coord in cand_coords` plus, *each iteration*,
   `np.cos(cand_coord.dec).value`, `cand_coord.ra.deg`,
   `cand_coord.dec.deg`.

**Bottleneck identified:** the per-candidate loop (item 4). Two costs:
- Iterating a SkyCoord array (`for cand_coord in cand_coords`) builds a
  fresh *scalar* SkyCoord object every candidate — expensive astropy
  object construction.
- The three attribute accesses per candidate each go through astropy's
  `Longitude`/`Latitude`/`Quantity` machinery.

Together this is ~100 us per candidate of pure astropy overhead that
has nothing to do with the actual math. I benchmarked it (N=50
candidates, 200x200 grid, best of 20):

| approach                              | time/call |
|---------------------------------------|-----------|
| iterate SkyCoord + attr per candidate | 8.75 ms   |
| pre-extract arrays once + index       | 3.12 ms   |
| (pure astropy attr access alone)      | 5.06 ms   |

So ~5.6 ms/call (2.8x on the theta-build slice) is removed simply by
pulling the candidate coordinates out of the SkyCoord *once* before the
loop. The saving is fixed per candidate, so it grows with the number of
candidates and dominates when grids are small.

Note: the genuinely heavy numeric work — building `theta`, `pw_Oi`
(`exp`/`log`), the `L_wx * p_wOi` product and `np.sum` over ngrid^2,
done N_cand times — is already pure numpy. That is O(N_cand * ngrid^2)
and is the floor for large grids; the astropy cleanup does not touch it.

**How to replace all astropy calculations:**
1. Before the loop, extract candidate coords to plain numpy arrays once:
   `cand_ra = cand_coords.ra.deg`, `cand_dec = cand_coords.dec.deg`,
   and precompute `cos_cand_dec = np.cos(np.radians(cand_dec))`. Inside
   the loop, index with `cand_ra[icand]` etc. — no scalar SkyCoord, no
   per-iteration astropy attribute access.
2. Extract the center once with numpy too: `center_ra = ...ra.deg`,
   `center_dec = ...dec.deg`, `cos_center_dec =
   np.cos(np.radians(center_dec))` (replaces the astropy
   `np.cos(...).value` on line 120).
3. Drop the equinox line (105). It is a no-op for the numpy paths
   (calc_LWx eellipse ignores equinox; healpix/wcs never used it) and
   the offsets here are computed flat-sky, not via spherical frames.

**On storing ra/dec in radians (the doc's hint):** modest, situational.
- The flat-sky `theta` (lines 129-130) uses *degree differences* scaled
  by 3600 to arcsec; in radians it would be `206265 * sqrt(...)`.
  Mathematically identical, no speed difference — not worth changing.
- The real redundancy is that `calc_LWx`'s eellipse branch now does
  `np.radians(ra)` and `np.radians(dec)` over the *whole* grid on every
  call. In `px_Oi_fixedgrid` that is a single call, so it is a one-time
  cost — minor here (it matters more for `px_Oi_local`, deferred). If we
  ever want to shave it, the clean move is to build the grid in radians
  up front and let `calc_LWx` accept radians (or cache the radian grid),
  rather than threading radians through the theta loop.

**Recommendations (no code written yet):**
1. Pre-extract candidate and center coordinates to numpy arrays once;
   index by `icand` in the loop. Biggest, safest win; output identical.
2. Remove the equinox-setting line.
3. (Optional, separate change) Vectorize the candidate loop with
   broadcasting — compute `theta`/`p_wOi` for all candidates at once as
   an (N, ngrid, ngrid) array. Fast for moderate sizes but memory grows
   as N_cand * ngrid^2 (e.g. 50 cands x 7200^2 x 8 B ~ 20 GB). Must keep
   the per-candidate loop for very large grids, or chunk it; flag the
   trade-off rather than vectorize unconditionally.
4. (Optional) `pw_Oi` already masks to `theta < max*phi`; `theta`
   itself is still computed over the full grid. Could restrict the
   sqrt/exp to the in-support region, but that complicates the code for
   limited gain — lower priority than items 1-2.

**Question for next prompt:** For item 1, confirm `cand_coords` is
always an array-valued SkyCoord (it is built that way in `path.py`) so
`cand_coords.ra.deg` returns an array — any single-candidate call sites
I should guard against? And do you want me to implement just items 1-2
(safe, identical output), or also pursue the vectorization in item 3
with a memory guard?

### 2026-06-09 (Implemented numpy px_Oi_fixedgrid items 1-2 + tests)

**Scope:** implemented only items 1 and 2 from the recommendations
above (pre-extract coords; drop equinox). No vectorization (item 3) or
`pw_Oi` masking (item 4), per the prompt.

**Code changed:** `astropath/bayesian.py:px_Oi_fixedgrid`.
- Item 2: removed `localiz['center_coord'].equinox = cand_coords[0]
  .equinox` (no-op for the numpy/flat-sky paths).
- Item 1 (center): extract `center_ra`, `center_dec` as plain floats
  once and precompute `cos_center_dec = np.cos(np.radians(center_dec))`,
  replacing the astropy `np.cos(...).value` in the grid construction.
- Item 1 (candidates): before the loop, pull
  `cand_ra = cand_coords.ra.deg`, `cand_dec = cand_coords.dec.deg`, and
  `cos_cand_dec = np.cos(np.radians(cand_dec))` as numpy arrays. The
  loop is now `for icand in range(cand_ra.size)` and indexes those
  arrays — no per-candidate scalar `SkyCoord` construction, no
  per-candidate astropy attribute access. Output is numerically
  identical (same flat-sky formula, just float vs Quantity extraction).

**numba consideration:** I kept the inner loop operating purely on
numpy arrays and scalar floats (`cand_ra[icand]`, `cos_cand_dec[icand]`,
`cand_ang_size[icand]`) with no astropy objects or dict access on the
hot path's coordinate math. That is exactly the shape `@njit` wants.
The remaining numba blocker is `pw_Oi`, which branches on the
`theta_prior` dict (string `PDF` key) — to jit the loop later that
PDF selection would need to be hoisted out (resolve the PDF + its
scalar norm once, pass plain floats in), but no restructuring was done
now since the prompt limited scope to items 1-2.

**Tests added** to `astropath/tests/test_speed_up.py`:
- `_px_Oi_fixedgrid_astropy(...)` — frozen verbatim copy of the original
  astropy loop/grid logic (returns the p(x|O_i) array), used as the
  regression + speed reference.
- `_make_candidates(...)` — deterministic helper (no RNG) building a
  candidate `SkyCoord` + angular sizes around a center.
- `test_px_Oi_fixedgrid_matches_astropy` — 10 candidates, asserts
  `np.allclose(rtol=1e-10, atol=1e-15)` vs the astropy reference.
- `test_px_Oi_fixedgrid_speed_up` — 50 candidates, 200x200 grid; times
  both, re-checks equality, prints the speed-up, asserts numpy is not
  slower.

**Measured result** (astro14 env, best-of-3, `pytest -s`):

| case                          | astropy  | numpy    | speed-up |
|-------------------------------|----------|----------|----------|
| 50 cand, 200x200 grid         | 17.19 ms | 10.23 ms | 1.7x     |

The ~1.7x reflects removing the fixed per-candidate astropy overhead;
the residual time is the genuine numpy compute (theta/pw_Oi/product/sum
over the grid, x50), which is unchanged and would be the next target
for numba or vectorization.

**Verification:** the 2 new tests pass (plus the 4 existing
test_speed_up tests); `test_path.py`, `test_bayesian.py`,
`test_localization.py` (7 tests) still pass — the fixed-grid posteriors
are unchanged downstream.

### 2026-06-09 (Profiling module for calc_LWx + px_Oi_fixedgrid)

**Code added:** `astropath/profiling.py`, modeled on
`calculations/step_size/Profiling.ipynb`. It reuses the notebook's
faux-FRB scenario: FRB at 21h44m25.255s -40d54m00.10s, a circular 5"
error ellipse, two candidate galaxies 1" away, `theta_prior =
dict(max=6, PDF='exp', scale=1)`, `box_hwidth=90"`.

Structure:
- `default_setup()` — builds the localiz / candidates / prior.
- `_build_grid()` — constructs the fixed ra/dec grid exactly as
  `px_Oi_fixedgrid` does (so `calc_LWx` is timed on a realistic grid).
- `_time_call()` — best-of-reps wall timer; reps scale down for big
  grids (1 rep > 1e6 px, 2 reps > 1e5, else 3) to keep the sweep quick.
- `run_profiling()` — sweeps `step_size` and returns a pandas table.
- `plot_results()` — log-log figure of time vs grid pixels.
- `main()` / `python -m astropath.profiling` — prints the table and
  writes `astropath/profiling_timing.png`.

**Sweep:** `step_size` in [0.5, 0.25, 0.1, 0.05, 0.025] arcsec which,
with box_hwidth=90", gives side lengths 360 -> 7200 px (the notebook's
max grid). Both `calc_LWx` and `px_Oi_fixedgrid` are timed at each size.

**Results** (astro14 env, printed to screen):

| step_size | ngrid | n_pixels   | calc_LWx (ms) | px_Oi_fixedgrid (ms) |
|-----------|-------|------------|---------------|----------------------|
| 0.500     | 360   | 129,600    | 8.37          | 10.24                |
| 0.250     | 720   | 518,400    | 40.16         | 52.89                |
| 0.100     | 1800  | 3,240,000  | 286.60        | 393.87               |
| 0.050     | 3600  | 12,960,000 | 1516.86       | 2337.06              |
| 0.025     | 7200  | 51,840,000 | 6230.52       | 9482.16              |

Both scale linearly with pixel count (slope ~1 on the log-log figure),
as expected for elementwise grid ops. `calc_LWx` is the bulk of
`px_Oi_fixedgrid`'s cost; the gap between the two curves is the
per-candidate theta/pw_Oi/product/sum work (here 2 candidates). Figure
saved to `astropath/profiling_timing.png`.

**Note:** with the prior astropy implementation the notebook recorded
the localization step at ~8.5 s for the 7200 grid (separation ~2.9 s +
PA ~3 s + coord setup/x,y/L_wx); the numpy `calc_LWx` now does the same
7200 grid in ~6.2 s end-to-end — consistent with the ~2x large-grid
speed-up measured earlier, and the astropy separation/PA cost is gone.