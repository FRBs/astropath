"""
Example usage of the assign_host module.

This script demonstrates how to:
1. Generate a population of simulated FRBs
2. Load or create a galaxy catalog
3. Assign FRBs to host galaxies by magnitude matching
4. Save and analyze the results
"""

import numpy as np
import pandas as pd
from astropath.simulations import generate_frbs, assign_frbs_to_hosts


def main():
    """Run a complete FRB-to-host assignment example."""

    print("=" * 70)
    print("Example: Assigning FRBs to Host Galaxies")
    print("=" * 70)

    # =========================================================================
    # Step 1: Generate FRB Population
    # =========================================================================
    print("\nStep 1: Generating FRB population")
    print("-" * 70)

    n_frbs = 500
    survey = 'CHIME'
    seed = 42

    print(f"Generating {n_frbs} FRBs for {survey} survey...")
    frbs = generate_frbs(n_frbs=n_frbs, survey=survey, seed=seed)

    print(f"  Generated {len(frbs)} FRBs")
    print(f"  DM range: {frbs['DM'].min():.1f} - {frbs['DM'].max():.1f} pc/cm³")
    print(f"  Redshift range: {frbs['z'].min():.3f} - {frbs['z'].max():.3f}")
    print(f"  Host M_r range: {frbs['M_r'].min():.2f} - {frbs['M_r'].max():.2f}")
    print(f"  Host m_r range: {frbs['m_r'].min():.2f} - {frbs['m_r'].max():.2f}")

    # Optional: Save FRBs to file
    frb_file = 'frbs_example.csv'
    frbs.to_csv(frb_file, index=False)
    print(f"  Saved to: {frb_file}")

    # =========================================================================
    # Step 2: Load or Create Galaxy Catalog
    # =========================================================================
    print("\nStep 2: Preparing galaxy catalog")
    print("-" * 70)

    # For this example, we'll create a mock catalog
    # In practice, you would load a real survey catalog (Pan-STARRS, DECaL, etc.)

    print("Creating mock galaxy catalog...")
    galaxies = create_mock_galaxy_catalog(n_galaxies=10000, seed=seed)

    print(f"  Created {len(galaxies)} galaxies")
    print(f"  RA range: {galaxies['ra'].min():.3f} - {galaxies['ra'].max():.3f} deg")
    print(f"  Dec range: {galaxies['dec'].min():.3f} - {galaxies['dec'].max():.3f} deg")
    print(f"  Magnitude range: {galaxies['mag_best'].min():.2f} - {galaxies['mag_best'].max():.2f}")

    # =========================================================================
    # Step 3: Assign FRBs to Host Galaxies
    # =========================================================================
    print("\nStep 3: Assigning FRBs to host galaxies")
    print("-" * 70)

    # Define localization error ellipse
    # For CHIME-like localization: ~0.5" x 0.3" ellipse
    localization = (0.5, 0.3, 45.)  # (a, b, PA) in arcsec, arcsec, degrees

    print(f"Localization error ellipse:")
    print(f"  Semi-major axis (a): {localization[0]}\"")
    print(f"  Semi-minor axis (b): {localization[1]}\"")
    print(f"  Position angle: {localization[2]}°")

    # Assign FRBs to hosts
    assignments = assign_frbs_to_hosts(
        frb_df=frbs,
        galaxy_catalog=galaxies,
        localization=localization,
        mag_range=(17., 28.),  # Only assign FRBs with hosts in this mag range
        scale=2.,  # Scale factor for galaxy half-light radius
        seed=seed
    )

    print(f"\nSuccessfully assigned {len(assignments)} FRBs to hosts")

    # =========================================================================
    # Step 4: Save Results
    # =========================================================================
    print("\nStep 4: Saving results")
    print("-" * 70)

    outfile = 'frb_host_assignments.csv'
    assignments.to_csv(outfile, index=False)
    print(f"  Saved assignments to: {outfile}")

    # =========================================================================
    # Step 5: Analyze Results
    # =========================================================================
    print("\nStep 5: Analyzing results")
    print("-" * 70)

    print("\nOffset Statistics:")
    print(f"  Galaxy offset:")
    print(f"    Mean: {assignments['gal_off'].mean():.3f}\"")
    print(f"    Median: {assignments['gal_off'].median():.3f}\"")
    print(f"    Std: {assignments['gal_off'].std():.3f}\"")

    print(f"  Localization offset:")
    print(f"    Mean: {assignments['loc_off'].mean():.3f}\"")
    print(f"    Median: {assignments['loc_off'].median():.3f}\"")
    print(f"    Std: {assignments['loc_off'].std():.3f}\"")

    print("\nHost Magnitude Statistics:")
    print(f"  Mean: {assignments['mag'].mean():.2f}")
    print(f"  Median: {assignments['mag'].median():.2f}")
    print(f"  Range: {assignments['mag'].min():.2f} - {assignments['mag'].max():.2f}")

    print("\n" + "=" * 70)
    print("Example completed successfully!")
    print("=" * 70)

    # Display sample results
    print("\nSample results (first 10 assignments):")
    print(assignments[['FRB_ID', 'ra', 'dec', 'gal_ID', 'mag', 'gal_off', 'loc_off']].head(10))

    return assignments


def create_mock_galaxy_catalog(n_galaxies=10000, seed=42):
    """
    Create a mock galaxy catalog for demonstration.

    In practice, you would load a real catalog from Pan-STARRS, DECaL, etc.

    Args:
        n_galaxies: Number of galaxies to generate
        seed: Random seed for reproducibility

    Returns:
        pd.DataFrame: Mock galaxy catalog
    """
    np.random.seed(seed)

    # Create galaxies in a 2x2 degree field
    ra = np.random.uniform(150.0, 152.0, n_galaxies)
    dec = np.random.uniform(2.0, 4.0, n_galaxies)

    # Magnitude distribution roughly matching deep survey depths
    # Using a distribution that peaks around m_r ~ 21-22
    mag_best = np.random.beta(a=2, b=5, size=n_galaxies) * 8 + 18
    mag_best = np.clip(mag_best, 18., 26.)

    # Half-light radii in arcsec (log-normal distribution)
    # Typical values: 0.1" - 3.0"
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


if __name__ == '__main__':
    assignments = main()

    print("\n" + "=" * 70)
    print("Next steps:")
    print("=" * 70)
    print("1. Load your real galaxy catalog (Pan-STARRS, DECaL, etc.)")
    print("2. Run PATH analysis on the assigned FRBs")
    print("3. Compare assigned hosts with PATH posteriors")
    print("4. Analyze association success rates and biases")
