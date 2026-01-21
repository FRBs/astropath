# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**astropath** is the Probabilistic Association of Transients to Hosts (PATH) - a Bayesian framework for associating astronomical transients (primarily Fast Radio Bursts/FRBs) with host galaxies. The package calculates posterior probabilities P(O_i|x) that a given galaxy candidate O_i is the true host of a transient at localized position x.

## Installation and Setup

```bash
# Install in development mode
pip install -e .

# Install with dependencies
pip install -r requirements.txt
```

## Running Tests

```bash
# Run all tests
python setup.py test

# Run specific test file
pytest astropath/tests/test_bayesian.py

# Run specific test function
pytest astropath/tests/test_bayesian.py::test_pw_Oi

# Run with verbose output
pytest -v astropath/tests/
```

## Building Documentation

```bash
cd docs/
make html
# Output will be in docs/_build/html/
```

## Code Architecture

### Core Workflow

The PATH analysis follows this conceptual flow:

1. **Localization** ([localization.py](astropath/localization.py)): Define the transient's probabilistic sky position
   - Supports error ellipses (`eellipse`), WCS maps (`wcs`), or HEALPix maps (`healpix`)
   - Key function: `calc_LWx()` computes L(w-x) for localization grid

2. **Candidate Selection** ([candidates.py](astropath/candidates.py)): Define potential host galaxies
   - Requires: RA, Dec, angular size (`ang_size`), and optionally apparent magnitude
   - Validation via `vet_candidates()`

3. **Priors** ([priors.py](astropath/priors.py)): Set prior probabilities
   - **Candidate priors** P(O_i): Methods include `inverse`, `inverse_ang`, `inverse_ang2`, `identical`, or user-defined
   - **Offset priors** p(θ|O_i): Distribution of transient offsets from galaxy centers
     - PDFs: `exp` (exponential), `core`, or `uniform`
     - Parameters: `max` (cutoff), `scale` (multiplicative factor)

4. **Bayesian Calculation** ([bayesian.py](astropath/bayesian.py)): Compute posteriors
   - `pw_Oi()`: Calculates p(w|O_i), the offset distribution for a galaxy
   - `px_Oi_fixedgrid()`: Computes p(x|O_i) by convolving localization with offset distribution
   - Final posterior: P(O_i|x) ∝ P(O_i) × p(x|O_i)

5. **PATH Class** ([path.py](astropath/path.py)): High-level interface orchestrating the workflow
   - Methods: `init_candidates()`, `init_cand_prior()`, `init_theta_prior()`, `init_localization()`, `calc_priors()`, `run_individual()`, `P_Oi_of_x()`

### Key Modules

- **[chance.py](astropath/chance.py)**: Galaxy number density calculations (Σ_m) using Driver et al. 2016 and Windhorst et al. magnitude counts
- **[montecarlo.py](astropath/montecarlo.py)**: Monte Carlo simulations for testing and validation with COSMOS data
- **[healpix.py](astropath/healpix.py)**: HEALPix localization handling, including conversion from error ellipses
- **[cosmos.py](astropath/cosmos.py)**: Interface to COSMOS survey data for simulations
- **[utils.py](astropath/utils.py)**: Common utilities for coordinate transformations and data handling

### Scripts

- **[bin/astropath_catalog](bin/astropath_catalog)**: Command-line tool to run PATH on a localization using public catalog data (Pan-STARRS or DECaL)
  - Implementation in [astropath/scripts/use_catalogs.py](astropath/scripts/use_catalogs.py)
  - Example: `astropath_catalog J081240.7+320809 0.5,0.3,45 --survey Pan-STARRS`

### Data Model Conventions

The package uses dictionaries with strict data models defined in the relevant modules:

- **`localization_dmodel`** ([localization.py](astropath/localization.py)): Keys include `type`, `center_coord`, `eellipse`, `healpix_data`, `wcs_data`, etc.
- **`theta_dmodel`** ([priors.py](astropath/priors.py)): Keys include `PDF`, `scale`, `max`
- **`cand_dmodel`** ([priors.py](astropath/priors.py)): Keys include `P_O_method`, `P_U`, `name`

Validation functions (e.g., `vet_theta_prior()`, `vet_cand_prior()`) enforce these models.

## Important Implementation Notes

- **Coordinate Systems**: All coordinates are ICRS unless specified. Error ellipse PA (`theta`) is East from North.
- **Units**: Angular sizes in arcsec, magnitudes assumed to be r-band (extinction-corrected), solid angles in steradians.
- **Grid Discretization**: Bayesian calculations use discretized grids with typical `step_size=0.1"` and `box_hwidth` defining the analysis region.
- **Normalization**: Prior probabilities are normalized via `renorm_priors()` to sum to (1 - P_U).
- **External Dependencies**: Requires `frb` package for survey utilities in catalog scripts.

## Common Patterns

### Setting up a PATH analysis

```python
from astropath import path

# Initialize
mypath = path.PATH()

# Define candidates
mypath.init_candidates(ra=ra_array, dec=dec_array, ang_size=size_array, mag=mag_array)

# Set priors
mypath.init_cand_prior(P_O_method='inverse', P_U=0.01)
mypath.init_theta_prior(PDF='exp', max=20, scale=0.5)

# Define localization
mypath.init_localization('eellipse', center_coord=coord, eellipse={'a': 1.0, 'b': 0.5, 'theta': 45})

# Run calculation
mypath.calc_priors()
P_Ox = mypath.calc_P_Ox(box_hwidth=10.)
```

### Adding a new offset PDF

1. Add new option to `theta_dmodel['PDF']['options']` in [priors.py](astropath/priors.py)
2. Implement the PDF in `pw_Oi()` in [bayesian.py](astropath/bayesian.py)
3. Ensure proper normalization (integrate to 1 over the support)
4. Add corresponding test in [astropath/tests/test_bayesian.py](astropath/tests/test_bayesian.py)

### Adding a new candidate prior method

1. Add method name to `cand_dmodel['P_O_method']['options']` in [priors.py](astropath/priors.py)
2. Implement calculation in `raw_prior_Oi()` in [priors.py](astropath/priors.py)
3. Add test in [astropath/tests/test_priors.py](astropath/tests/test_priors.py)

## Testing Guidelines

- Tests require pytest>=6.0.0
- Some tests are marked with `@remote_data` decorator and require `PATH_DATA` environment variable for extended test suites
- Test files mirror the module structure: `test_<module>.py` for each core module
- Normalization checks are critical for probability calculations
