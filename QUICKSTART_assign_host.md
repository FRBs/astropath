# Quick Start: assign_host Module

## Basic Usage

```python
from astropath.simulations import generate_frbs, assign_frbs_to_hosts

# 1. Generate FRBs
frbs = generate_frbs(1000, 'CHIME', seed=42)

# 2. Load galaxy catalog
import pandas as pd
galaxies = pd.read_csv('galaxy_catalog.csv')
# Required columns: ra, dec, mag_best, half_light, ID

# 3. Assign FRBs to hosts
assignments = assign_frbs_to_hosts(
    frbs, 
    galaxies,
    localization=(0.5, 0.3, 45.),  # a, b, PA in arcsec, deg
    seed=42
)

# 4. Save
assignments.to_csv('assignments.csv', index=False)
```

## Key Parameters

- `localization`: (a, b, PA) - error ellipse semi-major, semi-minor (arcsec), position angle (deg)
- `mag_range`: (min, max) - magnitude range for FRB filtering, default (17, 28)
- `scale`: float - scale factor for galaxy half-light radius, default 2.0
- `seed`: int - random seed for reproducibility

## Output Columns

- `ra`, `dec` - Observed coordinates (with localization error)
- `true_ra`, `true_dec` - True coordinates in galaxy
- `gal_ID` - Host galaxy ID
- `mag` - Galaxy magnitude
- `gal_off` - Offset from galaxy center (arcsec)
- `loc_off` - Localization error (arcsec)
- `FRB_ID` - Original FRB index

## Testing

```bash
# Run test
conda run -n astro python test_assign_host.py

# Run example
conda run -n astro python examples/assign_frbs_example.py
```

## Documentation

- Full docs: `astropath/simulations/README.md`
- Implementation details: `IMPLEMENTATION_SUMMARY.md`
- Code: `astropath/simulations/assign_host.py`
