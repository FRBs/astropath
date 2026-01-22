*************************
Quick Start: Host Assignment
*************************

This is a quick reference guide for the FRB-to-host assignment functionality
in ``astropath.simulations.assign_host``.

Basic Workflow
==============

The typical workflow is:

1. Generate FRBs with ``generate_frbs()``
2. Load or create a galaxy catalog
3. Assign FRBs to hosts with ``assign_frbs_to_hosts()``
4. Analyze and save results

Minimal Example
===============

::

    from astropath.simulations import generate_frbs, assign_frbs_to_hosts
    import pandas as pd

    # 1. Generate FRBs
    frbs = generate_frbs(1000, 'CHIME', seed=42)

    # 2. Load galaxy catalog
    galaxies = pd.read_csv('galaxy_catalog.csv')
    # Required columns: ra, dec, mag_best, half_light, ID

    # 3. Assign FRBs to hosts
    assignments = assign_frbs_to_hosts(
        frbs,
        galaxies,
        localization=(0.5, 0.3, 45.),  # a, b, PA in arcsec, deg
        seed=42
    )

    # 4. Save results
    assignments.to_csv('assignments.csv', index=False)

Key Parameters
==============

localization
------------

Error ellipse specification: ``(a, b, PA)``

* **a** -- Semi-major axis (arcsec)
* **b** -- Semi-minor axis (arcsec)
* **PA** -- Position angle (degrees, East of North)

Examples::

    # CHIME-like
    loc = (0.5, 0.3, 45.)

    # DSA-110 (high precision)
    loc = (0.1, 0.1, 0.)

    # ASKAP (larger)
    loc = (2.0, 1.0, 30.)

mag_range
---------

Magnitude range for FRB filtering: ``(min, max)``

Default: ``(17., 28.)``

Example::

    assignments = assign_frbs_to_hosts(
        frbs, galaxies, (0.5, 0.3, 45.),
        mag_range=(18., 26.),  # Only m_r in 18-26
        seed=42
    )

scale
-----

Scale factor for galaxy half-light radius: ``float``

Default: ``2.0``

* Smaller values = more concentrated near galaxy centers
* Larger values = more spread out

Example::

    assignments = assign_frbs_to_hosts(
        frbs, galaxies, (0.5, 0.3, 45.),
        scale=1.5,  # More concentrated
        seed=42
    )

seed
----

Random seed for reproducibility: ``int``

Example::

    assignments = assign_frbs_to_hosts(
        frbs, galaxies, (0.5, 0.3, 45.),
        seed=42  # Reproducible results
    )

Output Columns
==============

Coordinates
-----------

* ``ra``, ``dec`` -- Observed coordinates (with localization error)
* ``true_ra``, ``true_dec`` -- True coordinates in galaxy

Galaxy Properties
-----------------

* ``gal_ID`` -- Host galaxy ID
* ``mag`` -- Galaxy magnitude
* ``half_light`` -- Galaxy half-light radius (arcsec)

Offsets
-------

* ``gal_off`` -- Offset from galaxy center (arcsec)
* ``loc_off`` -- Localization error (arcsec)

Metadata
--------

* ``FRB_ID`` -- Original FRB index
* ``a``, ``b``, ``PA`` -- Localization parameters

Galaxy Catalog Format
=====================

Required columns in the galaxy catalog DataFrame:

========== ====== =====================================
Column     Type   Description
========== ====== =====================================
ra         float  Right ascension (degrees)
dec        float  Declination (degrees)
mag_best   float  Apparent r-band magnitude
half_light float  Half-light radius (arcseconds)
ID         int    Unique identifier
========== ====== =====================================

Example catalog structure::

    ra        dec       mag_best  half_light  ID
    150.234   2.456     21.5      0.45        0
    150.567   2.789     23.2      0.32        1
    151.123   3.001     19.8      0.78        2

File-Based Workflow
===================

Use ``assign_frbs_to_hosts_from_files()`` to read FRBs from a file::

    from astropath.simulations import assign_frbs_to_hosts_from_files

    assignments = assign_frbs_to_hosts_from_files(
        frb_file='frbs_chime.csv',
        galaxy_catalog=galaxies,
        localization=(0.5, 0.3, 45.),
        outfile='assignments.csv',
        seed=42
    )

This reads, assigns, and saves in one function call.

Common Patterns
===============

Multiple Surveys
----------------

::

    surveys = ['CHIME', 'DSA', 'ASKAP']
    localizations = {
        'CHIME': (0.5, 0.3, 45.),
        'DSA': (0.1, 0.1, 0.),
        'ASKAP': (2.0, 1.0, 30.)
    }

    for survey in surveys:
        frbs = generate_frbs(1000, survey, seed=42)
        assignments = assign_frbs_to_hosts(
            frbs, galaxies, localizations[survey], seed=42
        )
        assignments.to_csv(f'assignments_{survey}.csv', index=False)

Analyzing Results
-----------------

::

    import numpy as np

    # Offset statistics
    print(f"Mean galaxy offset: {assignments['gal_off'].mean():.3f}\"")
    print(f"Median galaxy offset: {assignments['gal_off'].median():.3f}\"")

    # Localization performance
    print(f"Mean localization error: {assignments['loc_off'].mean():.3f}\"")
    print(f"90th percentile: {np.percentile(assignments['loc_off'], 90):.3f}\"")

    # Host magnitude distribution
    print(f"Host magnitude range: {assignments['mag'].min():.1f} - "
          f"{assignments['mag'].max():.1f}")

Plotting
--------

::

    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    # Galaxy offsets
    axes[0].hist(assignments['gal_off'], bins=30, alpha=0.7)
    axes[0].set_xlabel('Galaxy Offset (arcsec)')
    axes[0].set_ylabel('Count')
    axes[0].set_title('FRB-Host Separations')

    # Localization errors
    axes[1].hist(assignments['loc_off'], bins=30, alpha=0.7)
    axes[1].set_xlabel('Localization Error (arcsec)')
    axes[1].set_ylabel('Count')
    axes[1].set_title('Localization Performance')

    plt.tight_layout()
    plt.savefig('offsets.png', dpi=150)

Testing
=======

Verify installation::

    # Basic import test
    from astropath.simulations import assign_frbs_to_hosts
    print("✓ Module imported successfully")

    # Quick functionality test
    frbs = generate_frbs(10, 'CHIME', seed=42)
    print(f"✓ Generated {len(frbs)} FRBs")

Run the example script::

    python examples/assign_frbs_example.py

Troubleshooting
===============

Common Issues
-------------

**"No FRBs remain after magnitude cut"**

* Solution: Adjust ``mag_range`` to match your FRB and galaxy distributions
* Check: ``print(frbs['m_r'].min(), frbs['m_r'].max())``

**"Missing required columns"**

* Solution: Ensure galaxy catalog has: ra, dec, mag_best, half_light, ID
* Check: ``print(galaxies.columns)``

**"FRBs could not be assigned"**

* Solution: Galaxy catalog may be too small or magnitude distribution mismatch
* Check: Compare FRB and galaxy magnitude distributions

Debug Mode
----------

Enable debug output for detailed information::

    assignments = assign_frbs_to_hosts(
        frbs, galaxies, (0.5, 0.3, 45.),
        debug=True,
        seed=42
    )

This prints iteration-by-iteration matching statistics.

Further Reading
===============

* :doc:`assign_host` -- Complete documentation
* :doc:`simulations` -- FRB generation
* :doc:`frb_example` -- Full PATH analysis workflow
