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

## Profiling

## Docs

## Prompts

1. Read this doc.  Proceed with the 1st item under Testing.

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
