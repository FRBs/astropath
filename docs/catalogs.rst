***************
Catalogs Module
***************

This doc describes the ``catalogs`` module which provides functionality
to query public astronomical surveys and prepare catalogs for PATH analysis.

Overview
========

The ``catalogs`` module provides functions for:

- Querying public surveys (Pan-STARRS, DECaL) for candidate host galaxies
- Cleaning catalogs to remove stars and problematic sources
- Adding NGC/IC galaxies from the OpenNGC database
- Applying dust corrections to magnitudes

This module requires the ``frb`` package for survey access. If the ``frb``
package is not installed, catalog querying will raise an ImportError with
instructions for installation.

Supported Surveys
=================

Pan-STARRS
----------

The Pan-STARRS (Panoramic Survey Telescope and Rapid Response System) survey
covers the entire sky north of declination -30°. When querying Pan-STARRS:

- Magnitude key: ``Pan-STARRS_r``
- Star/galaxy separation uses:
  - Kron radius cuts
  - PSF likelihood
  - PSF vs Kron magnitude comparison
  - Cross-match with Gaia DR3 to remove stars
  - Optional PS1-PSC (Point Source Catalog) classification

DECaL
-----

The DECaLS (Dark Energy Camera Legacy Survey) covers the extragalactic sky
visible from the northern hemisphere. When querying DECaL:

- Magnitude key: ``DECaL_r``
- Star/galaxy separation uses:
  - Type classification (removes PSF sources)
  - Magnitude cuts

Basic Usage
===========

Query a catalog for PATH analysis::

    from astropy.coordinates import SkyCoord
    from astropath import catalogs

    # Define transient position
    coord = SkyCoord('21h44m25.255s', '-40d54m00.10s', frame='icrs')

    # Query Pan-STARRS with 5 arcmin radius
    catalog, mag_key, stars = catalogs.query_catalog(
        'Pan-STARRS',
        coord,
        ssize=5.0,  # arcmin
        skip_NGC=False,
        dust_correct=True
    )

    print(f"Found {len(catalog)} candidates")
    print(f"Magnitude column: {mag_key}")
    print(f"Rejected {len(stars)} stars")

The returned catalog is an ``astropy.table.Table`` with the required
columns for PATH analysis:

- ``ra``: Right ascension (degrees)
- ``dec``: Declination (degrees)
- ``ang_size``: Angular size / half-light radius (arcsec)
- ``<mag_key>``: Magnitude in the survey band
- ``ID``: Source identifier

Integration with run Module
===========================

The ``catalogs`` module integrates with the ``run`` module for automatic
catalog querying. Set the ``survey`` key in the priors dictionary::

    from astropath import run

    # Build idict with survey for automatic querying
    idict = run.build_idict(
        ra=180.0,
        dec=-45.0,
        lparam={'a': 2.0, 'b': 1.0, 'theta': 0.},
        survey='Pan-STARRS',  # Will query Pan-STARRS automatically
        PU=0.1
    )

    # Run PATH - catalog is queried automatically
    result, P_Ux, Path, mag_key, cut_catalog, stars = run.run_on_dict(
        idict,
        verbose=True,
        skip_NGC=False,
        dust_correct=True
    )

This approach is convenient for quick analyses but requires network access
and the ``frb`` package.

Catalog Cleaning
================

The ``query_catalog()`` function applies several cleaning steps to remove
stars and problematic sources:

Pan-STARRS Cleaning
-------------------

1. **Non-detections**: Remove sources with r-band magnitude ≤ 0
2. **Size cut**: Remove sources with Kron radius ≤ 0
3. **Magnitude cut**: Remove sources brighter than r = 15 mag
4. **PSF likelihood**: Remove point-like sources (log(PSFLikelihood) ≥ -2)
5. **PSF-Kron comparison**: Remove sources where iPSFMag - iKronMag ≤ 0.05
6. **Gaia cross-match**: Remove sources matched to Gaia DR3 stars
7. **Faint exception**: Sources fainter than r = 20 mag are kept regardless
   (too faint for reliable star/galaxy classification)

DECaL Cleaning
--------------

1. **Magnitude cut**: Remove sources brighter than r = 14 mag or with NaN magnitudes
2. **Type cut**: Remove sources classified as PSF (point sources)
3. **Faint exception**: Sources fainter than r = 23 mag are kept regardless
4. **Zero size handling**: Sources with zero angular size are assigned 0.7 arcsec

Star/Galaxy Separation Parameters
=================================

For Pan-STARRS, you can provide custom star/galaxy separation parameters
via the ``star_galaxy_sep`` argument::

    from astropath import catalogs

    # Custom PS1-PSC cut (requires PS1-PSC data and function)
    star_galaxy_sep = {
        'PS1_PSC_cut': 0.5,      # More aggressive cut (default 0.83)
        'PS1_PSC_func': my_psc_function  # Custom function to get PSC values
    }

    catalog, mag_key, stars = catalogs.query_catalog(
        'Pan-STARRS',
        coord,
        5.0,
        star_galaxy_sep=star_galaxy_sep
    )

If ``PS1_PSC_func`` is not provided, the PS1-PSC cut is skipped.

NGC/IC Galaxies
===============

By default, ``query_catalog()`` adds bright NGC and IC galaxies within
60 arcmin of the transient position using the OpenNGC database::

    from astropy.coordinates import SkyCoord
    from astropy.table import Table
    from astropath import catalogs

    coord = SkyCoord(180.0, 45.0, unit='deg')
    catalog = Table()
    catalog['ra'] = [180.0]
    catalog['dec'] = [45.0]
    catalog['ang_size'] = [1.0]
    catalog['mag'] = [20.0]

    # Add NGC galaxies to existing catalog
    catalogs.add_NGC_galaxies(coord, catalog, 'mag', sep_arcmins=60)

The function modifies the catalog in place, adding rows for each NGC/IC
galaxy found. V-band magnitudes from OpenNGC are used as approximations
for r-band.

To skip NGC galaxy addition::

    catalog, mag_key, stars = catalogs.query_catalog(
        'Pan-STARRS',
        coord,
        5.0,
        skip_NGC=True  # Don't add NGC galaxies
    )

Dust Correction
===============

By default, magnitudes are corrected for Galactic dust extinction using
the ``frb`` package utilities::

    # Enable dust correction (default)
    catalog, mag_key, stars = catalogs.query_catalog(
        'Pan-STARRS',
        coord,
        5.0,
        dust_correct=True
    )

    # Disable dust correction
    catalog, mag_key, stars = catalogs.query_catalog(
        'Pan-STARRS',
        coord,
        5.0,
        dust_correct=False
    )

The dust correction uses the E(B-V) value at the first source position
and applies the appropriate extinction correction for the survey band.

Dependencies
============

Required
--------

- ``astropy``: For Table and SkyCoord operations
- ``astroquery``: For Vizier queries (Gaia cross-match)
- ``numpy``: For numerical operations

Optional
--------

- ``frb``: Required for survey queries (``frb.surveys``) and dust correction
  (``frb.galaxies.nebular``, ``frb.galaxies.photom``)
- ``pyongc``: Required for NGC/IC galaxy addition

If optional dependencies are not available, the affected functionality
will be skipped with a warning message.

Error Handling
==============

The module handles various error conditions gracefully:

- **Missing frb package**: Raises ImportError with installation instructions
- **Missing pyongc**: Prints warning and skips NGC galaxy addition
- **Failed Gaia query**: Prints warning and continues without Gaia filtering
- **Empty catalog**: Returns empty catalog with None for mag_key and stars
- **Unsupported survey**: Raises ValueError

Example::

    from astropath import catalogs

    try:
        catalog, mag_key, stars = catalogs.query_catalog(
            'SDSS',  # Not supported
            coord,
            5.0
        )
    except ValueError as e:
        print(f"Error: {e}")
        # Error: Not ready for this survey: SDSS

API Reference
=============

.. automodule:: astropath.catalogs
   :members:
   :undoc-members:
   :show-inheritance:
