#!/usr/bin/env python
"""
Quick test script for assign_host module.

This verifies the basic functionality of the FRB-to-host assignment code.
"""

import os
import numpy as np
import pandas as pd
from pathlib import Path
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


def test_real_galaxy_catalog():
    """
    Test assignment using real galaxy catalog if available.

    Uses the combined HSC/DECaLs/HECATE catalog if FRB_APATH is set.
    """
    print("\n" + "=" * 60)
    print("Testing with real galaxy catalog (if available)")
    print("=" * 60)

    # Check for FRB_APATH environment variable
    frb_apath = os.environ.get('FRB_APATH')
    if frb_apath is None:
        print("\n⚠ FRB_APATH not set, skipping real catalog test")
        return None

    catalog_path = Path(frb_apath) / 'combined_HSC_DECaLs_HECATE_galaxies_hecatecut.parquet'

    if not catalog_path.exists():
        print(f"\n⚠ Catalog not found at {catalog_path}")
        print("  Skipping real catalog test")
        return None

    print(f"\n✓ Found catalog at {catalog_path}")

    # Load the catalog
    print("\n1. Loading galaxy catalog...")
    galaxies = pd.read_parquet(catalog_path)
    print(f"   Loaded {len(galaxies)} galaxies")

    # Display catalog info
    print(f"\n   Catalog columns: {list(galaxies.columns)}")

    # Check for required columns
    required_cols = ['ra', 'dec', 'mag_best', 'half_light', 'ID']
    missing_cols = [col for col in required_cols if col not in galaxies.columns]

    if missing_cols:
        print(f"\n⚠ Missing required columns: {missing_cols}")
        print("  Available columns:", list(galaxies.columns))
        print("  Skipping test")
        return None

    print(f"   RA range: {galaxies['ra'].min():.2f} - {galaxies['ra'].max():.2f} deg")
    print(f"   Dec range: {galaxies['dec'].min():.2f} - {galaxies['dec'].max():.2f} deg")
    print(f"   Magnitude range: {galaxies['mag_best'].min():.2f} - {galaxies['mag_best'].max():.2f}")

    # Generate FRBs in the catalog footprint
    print("\n2. Generating FRBs...")
    frbs = generate_frbs(n_frbs=50, survey='CHIME', seed=42)
    print(f"   Generated {len(frbs)} FRBs")
    print(f"   Magnitude range: {frbs['m_r'].min():.2f} - {frbs['m_r'].max():.2f}")

    # Assign FRBs to hosts
    print("\n3. Assigning FRBs to hosts from real catalog...")
    localization = (0.5, 0.3, 45.)  # CHIME-like localization

    try:
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
        assert 'gal_ID' in assignments.columns, "Missing 'gal_ID' column"
        assert 'mag' in assignments.columns, "Missing 'mag' column"
        print("   ✓ All required columns present")

        # Check offsets
        print("\n5. Checking offset distributions...")
        print(f"   Galaxy offset: mean={assignments['gal_off'].mean():.3f}\", "
              f"median={assignments['gal_off'].median():.3f}\"")
        print(f"   Localization offset: mean={assignments['loc_off'].mean():.3f}\", "
              f"median={assignments['loc_off'].median():.3f}\"")
        print("   ✓ Offsets computed successfully")

        # Display sample results
        print("\n6. Sample assignments (first 5):")
        cols_to_show = ['ra', 'dec', 'gal_ID', 'mag', 'gal_off', 'loc_off']
        print(assignments[cols_to_show].head())

        print("\n" + "=" * 60)
        print("Real catalog test passed!")
        print("=" * 60)

        return assignments

    except Exception as e:
        print(f"\n✗ Error during assignment: {e}")
        import traceback
        traceback.print_exc()
        raise


if __name__ == '__main__':
    # Run basic test with mock catalog
    assignments = test_basic_assignment()
    print("\n✓ assign_host module is working correctly!")

    # Try real catalog test if available
    try:
        real_assignments = test_real_galaxy_catalog()
        if real_assignments is not None:
            print("\n✓ Real galaxy catalog test passed!")
    except Exception as e:
        print(f"\n✗ Real catalog test failed: {e}")
