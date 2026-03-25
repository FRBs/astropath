**********
Run Module
**********

This doc describes the ``run`` module which provides a high-level
interface for running PATH analysis on a single transient given
an input dictionary of parameters and a catalog of candidate host galaxies.

Overview
========

The ``run`` module provides three main functions:

- ``run_on_dict()``: Execute PATH analysis given configuration and catalog
- ``set_anly_sizes()``: Compute appropriate analysis sizes from localization
- ``build_idict()``: Convenience function to construct input dictionaries

This module is designed to be self-contained within astropath, with no
external dependencies beyond the standard astropath requirements.

Input Dictionary
================

The ``run_on_dict()`` function expects a dictionary with the following structure:

Required Keys
-------------

- **ra** (float): Right ascension in degrees
- **dec** (float): Declination in degrees
- **ssize** (float): Radius for probing the survey in arcmin
- **ltype** (str): Localization type (``'eellipse'`` or ``'healpix'``)
- **lparam** (dict): Localization parameters (see below)
- **priors** (dict): PATH prior configuration (see below)
- **max_box** (float): Maximum box size for PATH analysis in arcsec

Localization Parameters (lparam)
--------------------------------

For error ellipse (``ltype='eellipse'``)::

    lparam = {
        'a': 1.0,      # Semi-major axis in arcsec
        'b': 0.5,      # Semi-minor axis in arcsec
        'theta': 45.0  # Position angle in degrees, E from N
    }

For HEALPix (``ltype='healpix'``)::

    lparam = {
        'healpix_data': healpix_array,    # HEALPix probability map
        'healpix_nside': 512,             # NSIDE parameter
        'healpix_coord': 'C',             # Coordinate system ('C' for celestial)
        'healpix_ordering': 'RING'        # Pixel ordering ('RING' or 'NESTED')
    }

Prior Configuration (priors)
----------------------------

::

    priors = {
        'P_O_method': 'inverse',  # Method for candidate prior ['inverse', 'identical', etc.]
        'PU': 0.1,                # Prior probability of unseen host (0 to 1)
        'scale': 0.5,             # Scale factor for offset prior
        'theta_PDF': 'exp',       # PDF for offset prior ['exp', 'core', 'uniform']
        'theta_max': 6.0,         # Maximum offset for prior
        'survey': 'Pan-STARRS'    # Optional: survey for automatic catalog query
    }

The ``survey`` key is optional. If provided and no catalog is passed to
``run_on_dict()``, the catalog will be queried automatically using the
``catalogs`` module. See :doc:`catalogs` for details.

Candidate Catalog
=================

The catalog must be an ``astropy.table.Table`` with the following required columns:

- **ra**: Right ascension in degrees
- **dec**: Declination in degrees
- **ang_size**: Angular size (half-light radius) in arcsec
- **<mag_key>**: Magnitude column (name specified by ``mag_key`` parameter)

Optional columns:

- **ID**: Source identifier
- **z_phot_median**, **z_phot**, etc.: Photometric redshift columns (will be propagated to output)

Basic Usage
===========

Here is a complete example using the FRB example data::

    import os
    import pandas
    from importlib.resources import files as resource_files
    from astropy.table import Table
    from astropy.coordinates import SkyCoord

    from astropath import run

    # Load candidates
    cand_file = os.path.join(str(resource_files('astropath').joinpath('data')),
                             'frb_example', 'frb180924_candidates.csv')
    candidates_df = pandas.read_csv(cand_file, index_col=0)

    # Convert to astropy Table
    catalog = Table()
    catalog['ra'] = candidates_df.ra.values
    catalog['dec'] = candidates_df.dec.values
    catalog['ang_size'] = candidates_df.half_light.values
    catalog['mag'] = candidates_df.VLT_FORS2_g.values
    catalog['ID'] = candidates_df.index.values

    # FRB coordinates
    frb_coord = SkyCoord('21h44m25.255s -40d54m00.10s', frame='icrs')

    # Build input dictionary
    idict = run.build_idict(
        ra=frb_coord.ra.deg,
        dec=frb_coord.dec.deg,
        ltype='eellipse',
        lparam={'a': 0.1, 'b': 0.1, 'theta': 0.},
        P_O_method='inverse',
        PU=0.,
        scale=1.,
        theta_PDF='exp',
        theta_max=6.
    )

    # Run PATH
    result_candidates, P_Ux, Path, mag_key, cut_catalog, _ = run.run_on_dict(
        idict,
        verbose=True,
        catalog=catalog,
        mag_key='mag'
    )

The output ``result_candidates`` is a pandas DataFrame sorted by posterior
probability (``P_Ox``), with the most likely host first.

Using set_anly_sizes
====================

If you need to compute appropriate analysis sizes from localization
parameters, use ``set_anly_sizes()``::

    from astropath import run

    # For an error ellipse
    lparam = {'a': 5.0, 'b': 2.0, 'theta': 30.}
    ssize, max_box = run.set_anly_sizes('eellipse', lparam)

    print(f"Survey search radius: {ssize} arcmin")
    print(f"Analysis box size: {max_box} arcsec")

The function enforces a minimum search radius of 3 arcmin.

Using build_idict
=================

The ``build_idict()`` function simplifies creating input dictionaries::

    from astropath import run

    # Minimal usage - sizes are auto-computed
    idict = run.build_idict(
        ra=180.0,
        dec=-45.0,
        lparam={'a': 2.0, 'b': 1.0, 'theta': 0.}
    )

    # Full specification
    idict = run.build_idict(
        ra=180.0,
        dec=-45.0,
        ltype='eellipse',
        lparam={'a': 2.0, 'b': 1.0, 'theta': 0.},
        P_O_method='inverse',
        PU=0.1,
        scale=0.5,
        theta_PDF='exp',
        theta_max=6.0,
        ssize=5.0,      # Override auto-computed value
        max_box=300.0   # Override auto-computed value
    )

    # With automatic catalog querying
    idict = run.build_idict(
        ra=180.0,
        dec=-45.0,
        lparam={'a': 2.0, 'b': 1.0, 'theta': 0.},
        survey='Pan-STARRS'  # Will query Pan-STARRS when run_on_dict is called
    )

Automatic Catalog Querying
==========================

Instead of providing a catalog manually, you can let ``run_on_dict()``
query public surveys automatically. This requires the ``frb`` package.

Specify a survey in the input dictionary::

    from astropath import run

    # Build idict with survey
    idict = run.build_idict(
        ra=180.0,
        dec=45.0,
        lparam={'a': 5.0, 'b': 3.0, 'theta': 0.},
        survey='Pan-STARRS',
        PU=0.1
    )

    # Run PATH - catalog is queried automatically
    result, P_Ux, Path, mag_key, cut_catalog, stars = run.run_on_dict(
        idict,
        verbose=True,
        skip_NGC=False,      # Include NGC/IC galaxies
        dust_correct=True    # Apply dust correction
    )

    print(f"Using {mag_key} magnitudes")
    print(f"Rejected {len(stars) if stars else 0} stars")

Supported surveys:

- ``'Pan-STARRS'``: Pan-STARRS DR1 (dec > -30°)
- ``'DECaL'``: DECaLS DR10

See :doc:`catalogs` for detailed information about catalog querying,
cleaning, and star/galaxy separation.

Return Values
=============

The ``run_on_dict()`` function returns a tuple of 6 values:

1. **result_candidates** (pandas.DataFrame): Candidates with computed posteriors,
   sorted by ``P_Ox`` (descending). Contains columns:

   - ``ra``, ``dec``: Coordinates
   - ``ang_size``: Angular size
   - ``mag``: Magnitude
   - ``P_O``: Prior probability
   - ``P_Ox``: Posterior probability
   - ``P_Ux``: Probability of unseen host (same for all rows)
   - ``p_xO``: Likelihood p(x|O)
   - ``sep``: Separation from localization center (arcsec)
   - ``ID``: Source ID (if provided in input catalog)

2. **P_Ux** (float): Probability that the true host is unseen

3. **Path** (PATH): The PATH object used for analysis

4. **mag_key** (str): The magnitude column key used

5. **cut_catalog** (astropy.Table): Catalog cut to the analysis box

6. **stars** (astropy.Table or None): Table of sources rejected as stars
   during catalog querying. Only populated when using automatic catalog
   querying via the ``survey`` parameter; otherwise None.

Example with Unseen Prior
=========================

When you expect the true host might not be in the catalog, set ``PU > 0``::

    idict = run.build_idict(
        ra=frb_coord.ra.deg,
        dec=frb_coord.dec.deg,
        lparam={'a': 5.0, 'b': 5.0, 'theta': 0.},
        PU=0.2,  # 20% prior probability of unseen host
    )

    result, P_Ux, Path, _, _, _ = run.run_on_dict(
        idict, catalog=catalog, mag_key='mag'
    )

    # Check probability conservation
    total = result['P_Ox'].sum() + P_Ux
    assert abs(total - 1.0) < 1e-5  # Should sum to 1

API Reference
=============

.. automodule:: astropath.run
   :members:
   :undoc-members:
   :show-inheritance:
