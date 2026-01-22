#!/usr/bin/env python
"""
Quick test script for assign_host module.

This verifies the basic functionality of the FRB-to-host assignment code.
"""

import numpy as np
import pandas as pd
from astropath.simulations import generate_frbs, assign_frbs_to_hosts


def create_mock_galaxy_catalog(n_galaxies=5000, seed=42):
    """Create a mock galaxy catalog for testing."""
    np.random.seed(seed)

    # Create galaxies in a 1x1 degree field
    ra = np.random.uniform(150.0, 151.0, n_galaxies)
    dec = np.random.uniform(2.0, 3.0, n_galaxies)

    # Magnitude distribution roughly matching survey depths
    mag_best = np.random.uniform(18., 26., n_galaxies)

    # Half-light radii in arcsec (typical values)
    half_light = np.random.lognormal(mean=np.log(0.5), sigma=0.5, size=n_galaxies)
    half_light = np.clip(half_light, 0.1, 3.0)

    df = pd.DataFrame({
        'ra': ra,
        'dec': dec,
        'mag_best': mag_best,
        'half_light': half_light,
        'ID': np.arange(n_galaxies)
    })

    return df


def test_basic_assignment():
    """Test basic FRB-to-host assignment."""
    print("=" * 60)
    print("Testing basic FRB-to-host assignment")
    print("=" * 60)

    # Generate FRBs
    print("\n1. Generating FRBs...")
    frbs = generate_frbs(n_frbs=100, survey='CHIME', seed=42)
    print(f"   Generated {len(frbs)} FRBs")
    print(f"   Magnitude range: {frbs['m_r'].min():.2f} - {frbs['m_r'].max():.2f}")

    # Create galaxy catalog
    print("\n2. Creating mock galaxy catalog...")
    galaxies = create_mock_galaxy_catalog(n_galaxies=5000, seed=42)
    print(f"   Created {len(galaxies)} galaxies")
    print(f"   Magnitude range: {galaxies['mag_best'].min():.2f} - {galaxies['mag_best'].max():.2f}")

    # Assign FRBs to hosts
    print("\n3. Assigning FRBs to hosts...")
    localization = (0.5, 0.3, 45.)  # a=0.5", b=0.3", PA=45 deg

    assignments = assign_frbs_to_hosts(
        frbs,
        galaxies,
        localization=localization,
        mag_range=(17., 28.),
        seed=42
    )

    print(f"\n   Successfully assigned {len(assignments)} FRBs")

    # Validate results
    print("\n4. Validating results...")
    assert len(assignments) > 0, "No FRBs were assigned"
    assert 'ra' in assignments.columns, "Missing 'ra' column"
    assert 'dec' in assignments.columns, "Missing 'dec' column"
    assert 'true_ra' in assignments.columns, "Missing 'true_ra' column"
    assert 'gal_ID' in assignments.columns, "Missing 'gal_ID' column"
    assert 'mag' in assignments.columns, "Missing 'mag' column"
    assert 'loc_off' in assignments.columns, "Missing 'loc_off' column"

    # Check that localization parameters are stored
    assert np.all(assignments['a'] == localization[0])
    assert np.all(assignments['b'] == localization[1])
    assert np.all(assignments['PA'] == localization[2])

    print("   ✓ All required columns present")
    print("   ✓ Localization parameters stored correctly")

    # Check offsets are reasonable
    print("\n5. Checking offset distributions...")
    print(f"   Galaxy offset: mean={assignments['gal_off'].mean():.3f}\", "
          f"median={assignments['gal_off'].median():.3f}\"")
    print(f"   Localization offset: mean={assignments['loc_off'].mean():.3f}\", "
          f"median={assignments['loc_off'].median():.3f}\"")

    assert assignments['gal_off'].max() < 10., "Galaxy offsets suspiciously large"
    assert assignments['loc_off'].max() < 2., "Localization offsets suspiciously large"

    print("   ✓ Offsets are reasonable")

    # Check that observed coordinates differ from true coordinates
    ra_diff = np.abs(assignments['ra'] - assignments['true_ra'])
    dec_diff = np.abs(assignments['dec'] - assignments['true_dec'])
    total_diff = np.sqrt(ra_diff**2 + dec_diff**2) * 3600  # arcsec

    print(f"   Coordinate difference: mean={total_diff.mean():.3f}\", "
          f"max={total_diff.max():.3f}\"")
    assert np.all(total_diff > 0), "Some obs coords identical to true coords"
    print("   ✓ Localization error applied correctly")

    print("\n" + "=" * 60)
    print("All tests passed!")
    print("=" * 60)

    # Display sample results
    print("\nSample assignments (first 5):")
    print(assignments[['ra', 'dec', 'gal_ID', 'mag', 'gal_off', 'loc_off']].head())

    return assignments


if __name__ == '__main__':
    assignments = test_basic_assignment()
    print("\n✓ assign_host module is working correctly!")
