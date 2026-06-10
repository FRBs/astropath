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

## Prompts

1. Read this doc.  Proceed with the 1st item under Development/calc_LWx
2. Read this doc.  Proceed with the 2nd item under Development/calc_LWx
3. Read this doc.  Proceed with the 1st item under Development/px_Oi_fixedgrid


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