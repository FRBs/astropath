# Simulations Module

This subpackage provides tools for generating simulated FRB populations and assigning them to host galaxies for PATH analysis validation and testing.

## Modules

### `generate_frbs.py`

Generates simulated FRB populations with realistic properties (DM, z, M_r, m_r) for different surveys.

**Key function:** `generate_frbs(n_frbs, survey, dm_catalog=None, cosmo=None, seed=None)`

- Samples DM from survey-specific P(DM,z) grids
- Samples redshifts from P(z|DM) distributions
- Samples host galaxy absolute magnitudes from known FRB host distribution
- Computes apparent magnitudes from redshift and M_r

**Supported surveys:**
- CHIME
- DSA
- ASKAP / CRAFT
- Parkes
- FAST

### `assign_host.py`

Assigns simulated FRBs to host galaxies from a catalog by matching apparent magnitudes (m_r).

**Key function:** `assign_frbs_to_hosts(frb_df, galaxy_catalog, localization, ...)`

**Algorithm:**
1. Filters FRBs to reasonable magnitude range (default: 17-28)
2. Matches FRBs to galaxies by magnitude using a "fake coordinate" approach
   - Encodes magnitudes as declinations for efficient sky coordinate matching
   - Preferentially matches brighter FRBs to brighter galaxies
   - Each galaxy is assigned only once
3. Randomly places FRBs within their host galaxy (exponential offset distribution)
4. Applies localization error according to specified error ellipse

**Outputs:**
- Observed FRB coordinates (with localization error)
- True FRB coordinates (within galaxy)
- Galaxy ID and properties
- Offset statistics

**Based on:** `path_simulations.frbs.assign_chime_frbs_to_hosts()` and `frb.frb_surveys.mock.frbs_in_hosts()`

## Usage Example

```python
from astropath.simulations import generate_frbs, assign_frbs_to_hosts
import pandas as pd

# 1. Generate FRB population
frbs = generate_frbs(n_frbs=1000, survey='CHIME', seed=42)

# 2. Load galaxy catalog
# (must have columns: ra, dec, mag_best, half_light, ID)
galaxies = pd.read_csv('galaxy_catalog.csv')

# 3. Assign FRBs to hosts
# Localization: (semi-major, semi-minor, PA) in arcsec, arcsec, degrees
assignments = assign_frbs_to_hosts(
    frb_df=frbs,
    galaxy_catalog=galaxies,
    localization=(0.5, 0.3, 45.),
    mag_range=(17., 28.),
    seed=42
)

# 4. Save results
assignments.to_csv('frb_assignments.csv', index=False)
```

## Output Format

The `assign_frbs_to_hosts()` function returns a DataFrame with:

| Column | Description |
|--------|-------------|
| `ra` | Observed FRB RA (deg) - includes localization error |
| `dec` | Observed FRB Dec (deg) - includes localization error |
| `true_ra` | True FRB RA in galaxy (deg) |
| `true_dec` | True FRB Dec in galaxy (deg) |
| `gal_ID` | Assigned host galaxy ID |
| `gal_off` | Offset from galaxy center (arcsec) |
| `mag` | Galaxy apparent magnitude (m_r) |
| `half_light` | Galaxy half-light radius (arcsec) |
| `loc_off` | Localization error magnitude (arcsec) |
| `FRB_ID` | Original FRB index |
| `a` | Localization semi-major axis (arcsec) |
| `b` | Localization semi-minor axis (arcsec) |
| `PA` | Localization position angle (deg) |

## Galaxy Catalog Requirements

The galaxy catalog must be a pandas DataFrame with these columns:

- `ra`: Right ascension (degrees)
- `dec`: Declination (degrees)
- `mag_best`: Apparent r-band magnitude
- `half_light`: Half-light radius (arcsec)
- `ID`: Unique integer identifier

## Testing

Run the test suite:

```bash
python test_assign_host.py
```

Run the example:

```bash
python examples/assign_frbs_example.py
```

## Implementation Notes

### Magnitude Matching Algorithm

The assignment uses a clever "fake coordinate" approach:
1. Creates SkyCoord objects where declination = magnitude
2. Uses astropy's `match_coordinates_sky()` to match by "proximity" (= magnitude similarity)
3. Iteratively assigns FRBs, ensuring each galaxy is used only once
4. Handles edge cases where bright FRBs outnumber bright galaxies

This approach ensures:
- Bright FRBs → bright galaxies
- Faint FRBs → faint galaxies
- Realistic host associations for simulations

### Offset Distributions

**Galaxy offsets:** Truncated normal distribution scaled by galaxy half-light radius
- Generates offsets within 6σ of galaxy center
- Scale parameter controls concentration (default: 2.0)

**Localization error:** Truncated (3σ) normal along error ellipse axes
- Independent offsets along major (a) and minor (b) axes
- Total offset applied along position angle (PA)

### Coordinate System

All coordinates are ICRS. Position angles are measured East of North (IAU convention).

## Future Enhancements

Potential improvements:
- Support for asymmetric galaxy offset distributions
- Option for correlated localization errors
- Integration with COSMOS mock catalogs
- Built-in visualization tools
