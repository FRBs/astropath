***********************
Assigning FRBs to Hosts
***********************

This document describes the host galaxy assignment tools in *astropath*
for associating simulated FRBs with host galaxies.

Overview
========

The ``astropath.simulations.assign_host`` module provides tools for
assigning simulated Fast Radio Bursts (FRBs) to host galaxies from
a galaxy catalog. The assignment is based on matching apparent magnitudes
(m_r), ensuring that brighter FRBs are preferentially associated with
brighter galaxies.

This is useful for:

* Validating PATH analysis algorithms
* Testing host association methods
* Studying systematic biases in host identification
* Simulating realistic localization scenarios
* Planning survey strategies

The primary function is ``assign_frbs_to_hosts()`` which takes FRBs
from ``generate_frbs()`` and assigns them to galaxies, producing:

* **Observed coordinates**: FRB position including localization error
* **True coordinates**: Actual FRB position within the host galaxy
* **Galaxy properties**: Host galaxy ID, magnitude, size
* **Offset statistics**: Galaxy offset and localization error

Workflow
========

The typical workflow combines FRB generation with host assignment::

    from astropath.simulations import generate_frbs, assign_frbs_to_hosts
    import pandas as pd

    # Step 1: Generate FRB population
    frbs = generate_frbs(1000, 'CHIME', seed=42)

    # Step 2: Load galaxy catalog
    galaxies = pd.read_csv('galaxy_catalog.csv')

    # Step 3: Assign FRBs to hosts
    assignments = assign_frbs_to_hosts(
        frb_df=frbs,
        galaxy_catalog=galaxies,
        localization=(0.5, 0.3, 45.),  # a, b, PA in arcsec, deg
        seed=42
    )

    # Step 4: Save results
    assignments.to_csv('frb_assignments.csv', index=False)

Galaxy Catalog Requirements
============================

The galaxy catalog must be a pandas DataFrame with these required columns:

* **ra** (*float*) -- Right ascension in degrees (ICRS)
* **dec** (*float*) -- Declination in degrees (ICRS)
* **mag_best** (*float*) -- Apparent r-band magnitude
* **half_light** (*float*) -- Half-light radius in arcseconds
* **ID** (*int*) -- Unique integer identifier for each galaxy

Example catalog structure::

    ra        dec       mag_best  half_light  ID
    150.2341  2.4567    21.5      0.45        0
    150.5678  2.7890    23.2      0.32        1
    151.1234  3.0012    19.8      0.78        2
    ...

You can obtain suitable catalogs from:

* Pan-STARRS DR2
* DECaL/Legacy Survey
* SDSS
* Custom mock catalogs (e.g., COSMOS)

Basic Usage
===========

Simple Assignment
-----------------

The simplest usage assigns FRBs to a galaxy catalog with default
parameters::

    from astropath.simulations import generate_frbs, assign_frbs_to_hosts

    # Generate 100 CHIME FRBs
    frbs = generate_frbs(100, 'CHIME', seed=42)

    # Assign to hosts
    assignments = assign_frbs_to_hosts(
        frb_df=frbs,
        galaxy_catalog=galaxies,
        localization=(0.5, 0.3, 45.),
        seed=42
    )

    # View results
    print(assignments[['ra', 'dec', 'gal_ID', 'mag', 'gal_off', 'loc_off']].head())

The output will show:

* Observed FRB coordinates (``ra``, ``dec``)
* Assigned galaxy ID (``gal_ID``)
* Galaxy magnitude (``mag``)
* Offset from galaxy center (``gal_off``)
* Localization error (``loc_off``)

Localization Error Ellipse
---------------------------

The ``localization`` parameter specifies the error ellipse as a tuple
``(a, b, PA)`` where:

* **a** -- Semi-major axis in arcseconds
* **b** -- Semi-minor axis in arcseconds
* **PA** -- Position angle in degrees (East of North)

Examples for different surveys::

    # CHIME-like: ~0.5" x 0.3" ellipse
    loc_chime = (0.5, 0.3, 45.)

    # DSA-110: ~0.1" x 0.1" (nearly circular)
    loc_dsa = (0.1, 0.1, 0.)

    # ASKAP: larger ~1" x 2" ellipse
    loc_askap = (2.0, 1.0, 30.)

    assignments = assign_frbs_to_hosts(frbs, galaxies, loc_chime)

Magnitude Range Filtering
--------------------------

Control which FRBs are assigned by setting the magnitude range::

    # Only assign FRBs with hosts between 18-26 mag
    assignments = assign_frbs_to_hosts(
        frbs, galaxies,
        localization=(0.5, 0.3, 45.),
        mag_range=(18., 26.),
        seed=42
    )

FRBs outside this range are filtered out before assignment.

Galaxy Offset Distribution
---------------------------

The ``scale`` parameter controls how concentrated FRBs are near
galaxy centers::

    # More concentrated (closer to center)
    assignments_tight = assign_frbs_to_hosts(
        frbs, galaxies, (0.5, 0.3, 45.),
        scale=1.0,  # tighter distribution
        seed=42
    )

    # More spread out
    assignments_wide = assign_frbs_to_hosts(
        frbs, galaxies, (0.5, 0.3, 45.),
        scale=3.0,  # wider distribution
        seed=42
    )

The default ``scale=2.0`` produces realistic offsets based on observed
FRB host associations.

Output Format
=============

The ``assign_frbs_to_hosts()`` function returns a pandas DataFrame
with these columns:

Coordinates
-----------

* **ra** (*float*) -- Observed FRB RA in degrees (includes localization error)
* **dec** (*float*) -- Observed FRB Dec in degrees (includes localization error)
* **true_ra** (*float*) -- True FRB RA in the galaxy (degrees)
* **true_dec** (*float*) -- True FRB Dec in the galaxy (degrees)

Galaxy Properties
-----------------

* **gal_ID** (*int*) -- ID of the assigned host galaxy
* **mag** (*float*) -- Galaxy apparent r-band magnitude
* **half_light** (*float*) -- Galaxy half-light radius (arcsec)

Offsets
-------

* **gal_off** (*float*) -- Offset from galaxy center to true FRB position (arcsec)
* **loc_off** (*float*) -- Magnitude of localization error (arcsec)

Metadata
--------

* **FRB_ID** (*int*) -- Index of the FRB in the input DataFrame
* **a** (*float*) -- Localization semi-major axis (arcsec)
* **b** (*float*) -- Localization semi-minor axis (arcsec)
* **PA** (*float*) -- Localization position angle (degrees)

Example output::

    ra        dec       true_ra   true_dec  gal_ID  mag    gal_off  loc_off  FRB_ID
    150.4567  2.3456    150.4568  2.3457    1234    21.5   0.123    0.456    0
    150.7890  2.6789    150.7891  2.6790    5678    23.2   0.234    0.321    1
    ...

Algorithm Details
=================

Magnitude Matching
------------------

The core algorithm matches FRBs to galaxies by apparent magnitude using
a clever "fake coordinate" approach:

1. **Encode magnitudes**: Create SkyCoord objects where declination = magnitude
2. **Match coordinates**: Use astropy's ``match_coordinates_sky()`` to find nearest matches
3. **Iterative assignment**: Ensure each galaxy is used only once
4. **Handle edge cases**: Deal with situations where bright FRBs outnumber bright galaxies

This ensures realistic associations where brighter FRBs are preferentially
assigned to brighter galaxies, matching observed correlations.

Spatial Distributions
---------------------

**Galaxy Offsets**

FRB positions within galaxies are drawn from a truncated normal distribution:

* Centered on the galaxy center
* Width proportional to ``half_light / scale``
* Truncated at 6σ to avoid extreme outliers
* Random position angles

This produces realistic offset distributions comparable to observed
FRB-host separations.

**Localization Error**

Observed coordinates are offset from true positions following:

* Independent offsets along major (a) and minor (b) axes
* 3σ truncated normal distributions
* Combined via position angle rotation

The result mimics realistic telescope localization capabilities.

Advanced Usage
==============

Using a File for FRBs
---------------------

You can also read FRBs from a file using the convenience function::

    from astropath.simulations import assign_frbs_to_hosts_from_files

    # FRBs from a CSV file (output of generate_frbs)
    assignments = assign_frbs_to_hosts_from_files(
        frb_file='frbs_chime.csv',
        galaxy_catalog=galaxies,
        localization=(0.5, 0.3, 45.),
        outfile='assignments.csv',
        seed=42
    )

This reads the FRB file, assigns hosts, and saves the output in one call.

Custom Catalog Trimming
-----------------------

Control how much of the catalog edges are trimmed to maintain
valid PATH analysis regions::

    from astropy import units as u

    # Trim 2 arcmin from edges (default is 1 arcmin)
    assignments = assign_frbs_to_hosts(
        frbs, galaxies, (0.5, 0.3, 45.),
        trim_catalog=2*u.arcmin,
        seed=42
    )

This ensures FRBs don't fall too close to catalog boundaries.

Debug Mode
----------

Enable debug output to see detailed matching information::

    assignments = assign_frbs_to_hosts(
        frbs, galaxies, (0.5, 0.3, 45.),
        debug=True,
        seed=42
    )

This prints iteration-by-iteration matching statistics useful for
troubleshooting or understanding the algorithm.

Example: Complete Simulation
=============================

Here is a full example simulating CHIME FRBs and analyzing the
host associations::

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from astropath.simulations import generate_frbs, assign_frbs_to_hosts

    # Generate FRB population
    print("Generating FRBs...")
    frbs = generate_frbs(1000, 'CHIME', seed=42)

    # Load galaxy catalog (example: mock catalog)
    print("Loading galaxy catalog...")
    galaxies = pd.read_csv('panstarrs_catalog.csv')

    # Assign FRBs to hosts
    print("Assigning FRBs to hosts...")
    assignments = assign_frbs_to_hosts(
        frb_df=frbs,
        galaxy_catalog=galaxies,
        localization=(0.5, 0.3, 45.),
        mag_range=(17., 28.),
        seed=42
    )

    print(f"Successfully assigned {len(assignments)} FRBs")

    # Analyze results
    print("\\nOffset Statistics:")
    print(f"  Galaxy offset: {assignments['gal_off'].mean():.3f}\" ± "
          f"{assignments['gal_off'].std():.3f}\"")
    print(f"  Localization offset: {assignments['loc_off'].mean():.3f}\" ± "
          f"{assignments['loc_off'].std():.3f}\"")

    # Plot offset distributions
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    # Galaxy offsets
    axes[0].hist(assignments['gal_off'], bins=30, alpha=0.7, edgecolor='black')
    axes[0].set_xlabel('Galaxy Offset (arcsec)')
    axes[0].set_ylabel('Number of FRBs')
    axes[0].set_title('FRB Offsets from Galaxy Centers')
    axes[0].axvline(assignments['gal_off'].median(), color='red',
                    linestyle='--', label=f"Median: {assignments['gal_off'].median():.3f}\"")
    axes[0].legend()

    # Localization errors
    axes[1].hist(assignments['loc_off'], bins=30, alpha=0.7, edgecolor='black')
    axes[1].set_xlabel('Localization Error (arcsec)')
    axes[1].set_ylabel('Number of FRBs')
    axes[1].set_title('Localization Error Distribution')
    axes[1].axvline(assignments['loc_off'].median(), color='red',
                    linestyle='--', label=f"Median: {assignments['loc_off'].median():.3f}\"")
    axes[1].legend()

    plt.tight_layout()
    plt.savefig('frb_host_offsets.png', dpi=150)
    print("\\nPlot saved to: frb_host_offsets.png")

    # Save assignments
    assignments.to_csv('frb_host_assignments.csv', index=False)
    print("Assignments saved to: frb_host_assignments.csv")

API Reference
=============

assign_frbs_to_hosts()
----------------------

``assign_frbs_to_hosts(frb_df, galaxy_catalog, localization, mag_range=(17., 28.), scale=2., trim_catalog=1*u.arcmin, seed=None, debug=False)``

Assign FRBs to host galaxies based on apparent magnitude matching.

**Parameters:**

* **frb_df** (*pd.DataFrame*) -- FRB catalog with required column 'm_r'
* **galaxy_catalog** (*pd.DataFrame*) -- Galaxy catalog with columns: ra, dec, mag_best, half_light, ID
* **localization** (*tuple*) -- Error ellipse (a, b, PA) in arcsec, arcsec, degrees
* **mag_range** (*tuple, optional*) -- (min, max) magnitude range for FRB filtering (default: (17., 28.))
* **scale** (*float, optional*) -- Scale factor for galaxy half-light radius (default: 2.0)
* **trim_catalog** (*units.Quantity, optional*) -- Buffer to trim from catalog edges (default: 1 arcmin)
* **seed** (*int, optional*) -- Random seed for reproducibility
* **debug** (*bool, optional*) -- Enable debug output (default: False)

**Returns:**

* **pd.DataFrame** -- DataFrame with columns: ra, dec, true_ra, true_dec, gal_ID, mag, gal_off, loc_off, FRB_ID, a, b, PA

**Raises:**

* **ValueError** -- If required columns are missing from input DataFrames
* **RuntimeError** -- If FRBs cannot be assigned (catalog too small or magnitude mismatch)

assign_frbs_to_hosts_from_files()
----------------------------------

``assign_frbs_to_hosts_from_files(frb_file, galaxy_catalog, localization, outfile, **kwargs)``

Convenience function to assign FRBs from file and save results.

**Parameters:**

* **frb_file** (*str*) -- Path to CSV file containing FRB data
* **galaxy_catalog** (*pd.DataFrame*) -- Galaxy catalog DataFrame
* **localization** (*tuple*) -- Error ellipse (a, b, PA) in arcsec, arcsec, degrees
* **outfile** (*str*) -- Path to output CSV file
* **\\*\\*kwargs** -- Additional arguments passed to ``assign_frbs_to_hosts()``

**Returns:**

* **pd.DataFrame** -- Assignment results (also saved to outfile)

Implementation Notes
====================

Coordinate System
-----------------

* All coordinates are in ICRS
* Position angles are measured East of North (IAU convention)
* Units: degrees for coordinates, arcseconds for offsets

Random Number Generation
-------------------------

When a seed is specified, all random number generation is controlled
for full reproducibility. This includes:

* Magnitude-based galaxy selection
* Random placement within galaxies
* Localization error generation

Performance
-----------

Typical performance:

* 1000 FRBs with 10,000 galaxy catalog: ~1-2 seconds
* Memory efficient: operates on pandas DataFrames
* Scales approximately O(N log N) with number of FRBs

Dependencies
------------

Required packages:

* numpy, pandas (data handling)
* astropy (coordinates, units)
* scipy (used indirectly via coordinate matching)

The ``frb`` package is not required for ``assign_host`` itself,
but is needed if using ``generate_frbs()`` to create the FRB population.

See Also
========

* :doc:`simulations` -- Generating FRB populations
* :doc:`frb_example` -- Full PATH analysis example
* :doc:`localization` -- Localization methods in PATH
