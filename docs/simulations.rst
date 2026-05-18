***********
Simulations
***********

The ``astropath.simulations`` module provides a four-step pipeline for
testing and characterising PATH on synthetic FRB populations.

.. code-block:: text

   Step 1  generate_frbs()          →  FRB population (DMeg, z, M_r, m_r)
   Step 2  assign_frbs_to_hosts()   →  On-sky host assignments + localization
   Step 3  run_path.full()          →  PATH posteriors for every simulated FRB
   Step 4  utils.build_digest()     →  Merged results table ready for analysis

A command-line interface (:ref:`cli`) wraps all four steps into a single
invocation.  Interactive walkthroughs are available in the accompanying
notebooks (see :doc:`notebook_generate_frbs`, :doc:`notebook_assign_hosts`,
:doc:`notebook_run_path`).


.. _step1:

Step 1 — Generate an FRB Population
=====================================

``generate_frbs()`` samples a synthetic FRB population from a
survey-specific P(z, DM\ :sub:`EG`) grid provided by the ``frb`` package.
For each FRB it also draws a host-galaxy absolute magnitude from the known
FRB host luminosity function (``Lz_host_data.csv``), then converts to
apparent magnitude using the sampled redshift and a cosmological distance
modulus.

Supported surveys
-----------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Survey key
     - Description
   * - ``CHIME``
     - Canadian Hydrogen Intensity Mapping Experiment
   * - ``DSA``
     - Deep Synoptic Array (DSA-110)
   * - ``ASKAP``
     - Australian SKA Pathfinder — CRAFT class I & II
   * - ``CRAFT``
     - Alias for ``ASKAP``
   * - ``CRAFT_ICS_1300``
     - CRAFT Incoherent Sum at 1300 MHz
   * - ``CRAFT_ICS_892``
     - CRAFT Incoherent Sum at 892 MHz
   * - ``CRAFT_ICS_1632``
     - CRAFT Incoherent Sum at 1632 MHz
   * - ``Parkes``
     - Parkes Multibeam (class I & II)
   * - ``FAST``
     - Five-hundred-meter Aperture Spherical Telescope

API
---

``generate_frbs(n_frbs, survey, dm_catalog=None, cosmo=None, seed=None, dm_range=None)``

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Parameter
     - Default
     - Description
   * - *n_frbs*
     - —
     - Number of FRBs to generate *(int)*
   * - *survey*
     - —
     - Survey name; must be one of the keys above *(str)*
   * - *dm_catalog*
     - ``None``
     - 1-D array of observed extragalactic DMs.  When provided, DMs are
       drawn from a Gaussian KDE fitted to this array rather than
       directly from the P(DM, z) grid *(np.ndarray, optional)*
   * - *cosmo*
     - Planck18
     - Astropy cosmology for distance modulus calculations
       *(astropy.cosmology, optional)*
   * - *seed*
     - ``None``
     - Global random seed for reproducibility *(int, optional)*
   * - *dm_range*
     - ``None``
     - ``(min, max)`` DM range for the KDE evaluation when *dm_catalog*
       is provided *(tuple, optional)*

Output — ``pandas.DataFrame``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 15 10 75

   * - Column
     - Type
     - Description
   * - ``DMeg``
     - float
     - Extragalactic dispersion measure (pc cm\ :sup:`−3`)
   * - ``z``
     - float
     - Redshift
   * - ``M_r``
     - float
     - Host-galaxy absolute r-band magnitude
   * - ``m_r``
     - float
     - Host-galaxy apparent r-band magnitude

Usage examples
~~~~~~~~~~~~~~

Sample directly from the survey P(DM, z) grid::

    from astropath.simulations import generate_frbs

    frbs = generate_frbs(1000, 'CHIME', seed=42)

Constrain the DM distribution to an observed sample (e.g. CHIME
Catalog 1 bursts with S/N > 12)::

    from importlib.resources import files
    import pandas as pd
    import numpy as np

    fn = files('astropath.data') / 'frb_surveys' / 'chimefrbcat1.csv'
    df_dr1 = pd.read_csv(fn)
    dms_eg = np.nanmean([df_dr1['dm_exc_ne2001'].values,
                         df_dr1['dm_exc_ymw16'].values], axis=0)
    observed_dms = dms_eg[df_dr1['bonsai_snr'] > 12.]

    frbs = generate_frbs(1000, 'CHIME', dm_catalog=observed_dms, seed=42)


.. _step2:

Step 2 — Load a Galaxy Catalog and Assign FRBs to Hosts
=========================================================

Loading the host galaxy catalog
--------------------------------

``load_galaxy_catalog()`` reads a parquet catalog from the directory
given by the ``FRB_APATH`` environment variable::

    import os
    from astropath.simulations import load_galaxy_catalog

    os.environ['FRB_APATH'] = '/path/to/catalogs/'
    galaxies = load_galaxy_catalog(
        catalog_fn='combined_HSC_DECaLs_HECATE_galaxies_hecatecut.parquet'
    )

Two combined host catalogs are available from the
`Google Drive link <https://drive.google.com/drive/folders/1PKqh8tnDLbtqIuGeoPFEEh60ovs8Zjw8?usp=drive_link>`_
(~264 MB each).  Choose the one that matches the PATH catalog you will
use in Step 3:

* ``combined_HSC_PS1_HECATE_galaxies_hecatecut.parquet`` — pair with
  the Pan-STARRS PATH catalog
* ``combined_HSC_DECaLs_HECATE_galaxies_hecatecut.parquet`` — pair with
  the DECaLs PATH catalog

Both catalogs span HECATE (10 < *m*\ :sub:`r` < 14), Pan-STARRS or
DECaLs (14 < *m*\ :sub:`r` < 22), and HSC (22 < *m*\ :sub:`r` < 28).

The catalog must have the following columns:

.. list-table::
   :header-rows: 1
   :widths: 20 10 70

   * - Column
     - Type
     - Description
   * - ``ra``
     - float
     - Right ascension (degrees, ICRS)
   * - ``dec``
     - float
     - Declination (degrees, ICRS)
   * - ``mag``
     - float
     - Apparent r-band magnitude
   * - ``half_light``
     - float
     - Half-light radius (arcsec)
   * - ``ID``
     - int
     - Unique integer identifier

Assigning FRBs to host galaxies
--------------------------------

``assign_frbs_to_hosts()`` matches each FRB to a host galaxy by
apparent magnitude using a fake-coordinate iterative matching
algorithm, places the FRB within the host according to a chosen
intrinsic offset distribution, and offsets the position by a
localization error to produce the *observed* FRB coordinates.

API
~~~

``assign_frbs_to_hosts(frb_df, galaxy_catalog, localization, mag_range=None, offset_function='exponential', scale=0.5, trim_catalog=1*arcmin, seed=None, debug=False)``

.. list-table::
   :header-rows: 1
   :widths: 22 18 60

   * - Parameter
     - Default
     - Description
   * - *frb_df*
     - —
     - FRB DataFrame from ``generate_frbs()``; must contain column
       ``m_r`` *(pd.DataFrame)*
   * - *galaxy_catalog*
     - —
     - Galaxy catalog DataFrame; see column requirements above
       *(pd.DataFrame)*
   * - *localization*
     - —
     - Error ellipse as ``(a, b, PA)`` where *a* and *b* are the
       semi-major and semi-minor axes in arcsec and *PA* is the
       position angle in degrees East of North *(tuple)*
   * - *mag_range*
     - ``None``
     - ``(min, max)`` apparent magnitude range; FRBs outside this range
       are excluded.  ``None`` keeps all FRBs *(tuple, optional)*
   * - *offset_function*
     - ``'exponential'``
     - Intrinsic galactocentric offset distribution.
       One of ``'exponential'``, ``'uniform_1d'``, ``'uniform_2d'``
       *(str)*
   * - *scale*
     - ``0.5``
     - Scale parameter for the chosen offset distribution, in units of
       the host half-light radius *(float)*
   * - *trim_catalog*
     - ``1 arcmin``
     - Buffer trimmed from catalog edges to keep FRBs within the valid
       analysis region *(astropy.units.Quantity)*
   * - *seed*
     - ``None``
     - Random seed *(int, optional)*
   * - *debug*
     - ``False``
     - Print iteration-by-iteration matching diagnostics *(bool)*

Intrinsic offset distributions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - ``offset_function``
     - Distribution
   * - ``'exponential'``
     - :math:`p(r) \propto r\,\exp(-r/s)` — Gamma(2, scale) in units of
       the host half-light radius.  Matches the PATH exponential prior.
   * - ``'uniform_1d'``
     - Uniform in *r*, truncated at *scale* × half-light radius.
   * - ``'uniform_2d'``
     - Uniform per unit solid angle over a disk of radius *scale* ×
       half-light radius.  Matches the PATH uniform prior.

Output — ``pandas.DataFrame``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 18 10 72

   * - Column
     - Type
     - Description
   * - ``ra``
     - float
     - Observed FRB RA (degrees) — includes localization error
   * - ``dec``
     - float
     - Observed FRB Dec (degrees) — includes localization error
   * - ``true_ra``
     - float
     - True FRB RA inside the host galaxy (degrees)
   * - ``true_dec``
     - float
     - True FRB Dec inside the host galaxy (degrees)
   * - ``gal_ID``
     - int
     - ID of the assigned host galaxy in the catalog
   * - ``gal_off``
     - float
     - Angular offset from host-galaxy center to true FRB position (arcsec)
   * - ``mag``
     - float
     - Apparent r-band magnitude of the host galaxy
   * - ``half_light``
     - float
     - Half-light radius of the host galaxy (arcsec)
   * - ``loc_off``
     - float
     - Magnitude of the applied localization error (arcsec)
   * - ``FRB_ID``
     - int
     - Index of the FRB in the input ``frb_df``
   * - ``a``
     - float
     - Localization semi-major axis used (arcsec)
   * - ``b``
     - float
     - Localization semi-minor axis used (arcsec)
   * - ``PA``
     - float
     - Localization position angle used (degrees)

Usage example::

    from astropath.simulations import assign_frbs_to_hosts
    from astropy import units

    assignments = assign_frbs_to_hosts(
        frb_df          = frbs,
        galaxy_catalog  = galaxies,
        localization    = (25., 2., 12.),   # a, b, PA  [arcsec, arcsec, deg]
        offset_function = 'exponential',
        scale           = 0.5,
        trim_catalog    = 60 * units.arcmin,
        seed            = 42,
    )


.. _step3:

Step 3 — Run PATH
==================

``run_path.full()`` runs PATH on every simulated FRB scenario.  The
computation is embarrassingly parallel and is distributed across CPUs
with Python's ``multiprocessing`` module.  As a reference point: 5,000
FRBs took ~50 minutes on an Apple M3 MacBook Pro (16 GB RAM, ``ncpu = 6``).

Before calling this function, load the pre-queried PATH galaxy catalog.
Two options are available from the
`Google Drive link <https://drive.google.com/drive/folders/1PKqh8tnDLbtqIuGeoPFEEh60ovs8Zjw8?usp=drive_link>`_:

* ``catalog_dudxmmlss_hecate_Pan-STARRS.parquet`` (900 MB)
* ``catalog_dudxmmlss_hecate_DECaL.parquet`` (6.1 GB)

Choose the catalog that matches the host-galaxy catalog used in Step 2,
save it to ``$FRB_APATH``, and load it::

    path_catalog = load_galaxy_catalog(
        catalog_fn='catalog_dudxmmlss_hecate_DECaL.parquet'
    )

This catalog must have the following columns in addition to ``ra``,
``dec``, and ``ID``:

.. list-table::
   :header-rows: 1
   :widths: 20 10 70

   * - Column
     - Type
     - Description
   * - ``ang_size``
     - float
     - Angular size / half-light radius (arcsec)
   * - ``mag``
     - float
     - Apparent r-band magnitude

API
---

``run_path.full(frbs, catalog, prior_dict, multi=True, ncpu=4, debug=False)``

.. list-table::
   :header-rows: 1
   :widths: 20 12 68

   * - Parameter
     - Default
     - Description
   * - *frbs*
     - —
     - Host-assignment DataFrame from ``assign_frbs_to_hosts()``
       *(pd.DataFrame)*
   * - *catalog*
     - —
     - PATH galaxy catalog *(pd.DataFrame)*
   * - *prior_dict*
     - —
     - Dictionary of PATH prior settings; see below *(dict)*
   * - *multi*
     - ``True``
     - Use multiprocessing *(bool)*
   * - *ncpu*
     - ``4``
     - Number of parallel worker processes *(int)*
   * - *debug*
     - ``False``
     - Process only the first 100 FRBs *(bool)*

Prior dictionary
~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 22 20 58

   * - Key
     - Values
     - Description
   * - ``'P_O_method'``
     - ``'inverse'``, ``'identical'``
     - Magnitude prior :math:`P(O_i)`.  ``'inverse'`` weights by inverse
       flux; ``'identical'`` gives equal weight to all candidates.
   * - ``'PU'``
     - float in (0, 1)
     - Prior probability that the true host is unseen in the catalog.
       See Section 4.1.1 of Andersen+26 for guidance on estimating this
       value.
   * - ``'theta_PDF'``
     - ``'exp'``, ``'core'``, ``'uniform'``
     - Shape of the galactocentric offset prior
       :math:`p(\omega | O_i)`.  See :doc:`offset_function`.
   * - ``'scale'``
     - float
     - Multiplicative scale applied to the host half-light radius
       :math:`\phi` in the offset prior model.
   * - ``'theta_max'``
     - float
     - Maximum offset cutoff in units of :math:`\phi`; the prior is
       zero for :math:`\theta/\phi > \theta_\mathrm{max}`.

Output — ``pandas.DataFrame``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each row is the best PATH candidate for one simulated FRB.  Key columns:

.. list-table::
   :header-rows: 1
   :widths: 15 10 75

   * - Column
     - Type
     - Description
   * - ``ra``
     - float
     - RA of the candidate galaxy (degrees)
   * - ``dec``
     - float
     - Dec of the candidate galaxy (degrees)
   * - ``mag``
     - float
     - Apparent magnitude of the candidate galaxy
   * - ``ang_size``
     - float
     - Angular size of the candidate galaxy (arcsec)
   * - ``P_O``
     - float
     - PATH prior :math:`P(O_i)`
   * - ``p_xO``
     - float
     - PATH likelihood :math:`p(x | O_i)`
   * - ``P_Ox``
     - float
     - PATH posterior :math:`P(O_i | x)`
   * - ``P_Ux``
     - float
     - PATH unseen-host posterior :math:`P(U | x)`
   * - ``iFRB``
     - int
     - Index of the FRB in the ``frbs`` DataFrame

Usage example::

    from astropath.simulations import run_path

    prior_dict = {
        'P_O_method': 'inverse',
        'PU':         0.15,
        'theta_PDF':  'exp',
        'scale':      0.5,
        'theta_max':  6.0,
    }

    final_sims = run_path.full(
        frbs       = assignments,
        catalog    = path_catalog,
        prior_dict = prior_dict,
        multi      = True,
        ncpu       = 6,
    )


.. _step4:

Step 4 — Build the Digest
==========================

``utils.build_digest()`` merges the outputs of all three preceding steps
into a single, analysis-ready DataFrame and optionally writes it to a
parquet file.

API
---

``utils.build_digest(raw_sim_results, frbs, hosts, combined_catalog, output_fn=None, thresh_cross_match=2.)``

.. list-table::
   :header-rows: 1
   :widths: 25 12 63

   * - Parameter
     - Default
     - Description
   * - *raw_sim_results*
     - —
     - Output of ``run_path.full()`` *(pd.DataFrame)*
   * - *frbs*
     - —
     - Output of ``generate_frbs()`` *(pd.DataFrame)*
   * - *hosts*
     - —
     - Output of ``assign_frbs_to_hosts()`` *(pd.DataFrame)*
   * - *combined_catalog*
     - —
     - Host-galaxy catalog passed to ``assign_frbs_to_hosts()``
       *(pd.DataFrame)*
   * - *output_fn*
     - ``None``
     - If provided, the digest is written to this parquet path *(str)*
   * - *thresh_cross_match*
     - ``2.0``
     - A candidate is "correct" if its angular separation from the true
       host is within this multiple of the larger of the two galaxies'
       half-light radii *(float)*

Digest column reference
-----------------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Column
     - Description
   * - ``ra_loc``
     - Observed FRB RA — includes localization error (degrees)
   * - ``dec_loc``
     - Observed FRB Dec — includes localization error (degrees)
   * - ``true_ra``
     - True FRB RA inside the host galaxy (degrees)
   * - ``true_dec``
     - True FRB Dec inside the host galaxy (degrees)
   * - ``host_ID``
     - ID of the assigned host galaxy in the host catalog
   * - ``gal_off``
     - Offset from galaxy center to true FRB position (arcsec)
   * - ``mag_host``
     - Apparent r-band magnitude of the true host
   * - ``ang_size_host``
     - Angular size of the true host (arcsec)
   * - ``loc_off``
     - Localization error offset (arcsec)
   * - ``FRB_ID``
     - FRB index from ``generate_frbs()``
   * - ``a``
     - Localization semi-major axis (arcsec)
   * - ``b``
     - Localization semi-minor axis (arcsec)
   * - ``PA``
     - Localization position angle (degrees)
   * - ``ra_host``
     - RA of the host-galaxy center (degrees)
   * - ``dec_host``
     - Dec of the host-galaxy center (degrees)
   * - ``sep_best_host_arcsec``
     - Separation between the best candidate and the true host (arcsec)
   * - ``sep_host_loc_arcsec``
     - Separation between the true host center and the localization (arcsec)
   * - ``sep_best_loc_arcsec``
     - Separation between the best candidate center and the localization (arcsec)
   * - ``sep_host_loc_norm``
     - ``sep_host_loc_arcsec`` normalised by ``ang_size_host``
   * - ``sep_best_loc_norm``
     - ``sep_best_loc_arcsec`` normalised by ``ang_size_cand``
   * - ``z_host``
     - Simulated FRB redshift
   * - ``dmex_host``
     - Simulated FRB extragalactic DM (pc cm\ :sup:`−3`)
   * - ``frb_mr``
     - Simulated FRB host apparent r-band magnitude
   * - ``frb_Mr``
     - Simulated FRB host absolute r-band magnitude
   * - ``ra_cand``
     - RA of the best PATH candidate (degrees)
   * - ``dec_cand``
     - Dec of the best PATH candidate (degrees)
   * - ``mag_cand``
     - Apparent magnitude of the best candidate
   * - ``ang_size_cand``
     - Angular size of the best candidate (arcsec)
   * - ``cand_ID``
     - ID of the best candidate in the PATH catalog
   * - ``P_O``
     - PATH prior :math:`P(O_i)` for the best candidate
   * - ``p_xO``
     - PATH likelihood :math:`p(x | O_i)` for the best candidate
   * - ``P_Ox``
     - PATH posterior :math:`P(O_i | x)` for the best candidate
   * - ``P_Ux``
     - PATH unseen-host posterior :math:`P(U | x)`
   * - ``correct_association``
     - ``True`` if the best candidate spatially matches the true host
       (within *thresh_cross_match* × half-light radius)

Usage example::

    from astropath.simulations import utils as sim_utils

    digest = sim_utils.build_digest(
        raw_sim_results  = final_sims,
        frbs             = frbs,
        hosts            = assignments,
        combined_catalog = galaxies,
        output_fn        = 'digest.parquet',
    )

    # Quick performance summary
    thresh = 0.9
    n = len(digest)
    correct   = (digest['correct_association'] & (digest['P_Ox'] > thresh)).sum()
    incorrect = (~digest['correct_association'] & (digest['P_Ox'] > thresh)).sum()
    noassoc   = (digest['P_Ox'] < thresh).sum()
    print(f"Correct:     {correct}  ({100*correct/n:.1f}%)")
    print(f"Incorrect:   {incorrect}  ({100*incorrect/n:.1f}%)")
    print(f"Non-assoc.:  {noassoc}  ({100*noassoc/n:.1f}%)")


.. _cli:

Command-Line Interface
======================

``run_path_simulation_cli.py`` wraps the entire four-step pipeline into a
single command.  It enforces valid choices and types for every parameter,
generates output filenames automatically from the simulation settings, and
falls back gracefully to mock catalogs when real catalog files are not
present.

Usage
-----

.. code-block:: bash

   python run_path_simulation_cli.py \
       --survey          CHIME         \
       --n-frbs          1000          \
       --loc-a           25            \
       --loc-b           2             \
       --loc-pa          12            \
       --offset-dist     exponential   \
       --offset-scale    0.5           \
       --mag-prior       inverse       \
       --unseen-prior    0.15          \
       --offset-prior    exp           \
       --prior-scale     0.5           \
       --theta-max       6.0           \
       --ncpu            4             \
       --output-dir      ./results

Key arguments
~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Argument
     - Default
     - Description
   * - ``--survey``
     - required
     - Radio survey; see :ref:`step1` for valid choices
   * - ``--n-frbs``
     - 1000
     - Number of FRBs to simulate
   * - ``--loc-a / --loc-b / --loc-pa``
     - 25 / 2 / 12
     - Localization ellipse semi-axes (arcsec) and position angle (deg)
   * - ``--offset-dist``
     - ``exponential``
     - Intrinsic offset distribution (``exponential``, ``uniform_1d``,
       ``uniform_2d``)
   * - ``--offset-scale``
     - 0.5
     - Scale of the intrinsic offset distribution
   * - ``--mag-prior``
     - ``inverse``
     - PATH magnitude prior (``inverse``, ``identical``)
   * - ``--unseen-prior``
     - 0.15
     - PATH unseen-host prior P(U); must be in (0, 1)
   * - ``--offset-prior``
     - ``exp``
     - PATH offset prior shape (``exp``, ``core``, ``uniform``)
   * - ``--prior-scale``
     - 0.5
     - Scale parameter for the PATH offset prior
   * - ``--theta-max``
     - 6.0
     - Maximum offset cutoff for the PATH prior (units of :math:`\phi`)
   * - ``--ncpu``
     - 4
     - Number of parallel CPUs for PATH
   * - ``--host-catalog``
     - (see note)
     - Filename of the wide-field host catalog (looked up in
       ``$FRB_APATH``)
   * - ``--path-catalog``
     - (see note)
     - Filename of the pre-queried PATH galaxy catalog (looked up in
       ``$FRB_APATH``)
   * - ``--output-dir``
     - ``./``
     - Output directory
   * - ``--full-output``
     - auto
     - Parquet path for the raw PATH results
   * - ``--digest-output``
     - auto
     - Parquet path for the digest
   * - ``--seed``
     - ``None``
     - Global random seed
   * - ``--no-multi``
     - off
     - Disable multiprocessing (serial execution)
   * - ``--debug``
     - off
     - Process only 100 FRBs; useful for testing

Default catalog filenames (looked up in ``$FRB_APATH``):

* Host catalog: ``combined_HSC_DECaLs_HECATE_galaxies_hecatecut.parquet``
* PATH catalog: ``catalog_dudxmmlss_hecate_DECaL.parquet``

Auto-generated output filenames follow the pattern::

    digest__<loc>_<offset_dist>_<offset_prior>_<mag_prior>_pu<PU>_<catalog_stem>.parquet

so runs with different settings never overwrite each other.

Run ``python run_path_simulation_cli.py --help`` for the complete
argument list.


.. _andersen26:

Reproducing Andersen+26 Figures
================================

The notebook ``docs/nb/Reproducing_Andersen+26_Figures.ipynb`` contains
code to reproduce Figures 3–16 of Andersen+26.  To run it, download the
pre-computed simulation digests from the
`Google Drive link <https://drive.google.com/drive/folders/1Rv6koYfKJ7yQ5xYV666gIaMQ2wX0A6Tg?usp=drive_link>`_
and save them to ``$FRB_APATH``.  Digest files are named::

    digest_<loc>_<offset_dist>_<offset_prior>_<mag_prior>_<PU>_<catalog>.parquet

For example, ``digest_ck_exp_exp_inverse_0.15_DECaL.parquet`` corresponds
to the CHIME-KKO localization, exponential intrinsic offset distribution,
exponential PATH offset prior, inverse magnitude prior, P(U) = 0.15, run
with the DECaLs catalog (the first entry in Table 11 of Andersen+26).

Figure 7 (estimation of P(U)) additionally requires the pre-queried PATH
catalogs; see :ref:`step3` for download instructions.


Complete Example
================

The following reproduces a typical CHIME simulation end-to-end::

    import os
    import numpy as np
    import pandas as pd
    from importlib.resources import files
    from astropy import units
    from astropath.simulations import generate_frbs, assign_frbs_to_hosts
    from astropath.simulations import load_galaxy_catalog, run_path
    from astropath.simulations import utils as sim_utils

    os.environ['FRB_APATH'] = '/path/to/catalogs/'

    # Step 1 — generate FRBs
    fn = files('astropath.data') / 'frb_surveys' / 'chimefrbcat1.csv'
    df_dr1 = pd.read_csv(fn)
    dms_eg = np.nanmean([df_dr1['dm_exc_ne2001'].values,
                         df_dr1['dm_exc_ymw16'].values], axis=0)
    observed_dms = dms_eg[df_dr1['bonsai_snr'] > 12.]
    frbs = generate_frbs(1000, 'CHIME', dm_catalog=observed_dms, seed=42)

    # Step 2 — assign to hosts
    galaxies    = load_galaxy_catalog(
                      'combined_HSC_DECaLs_HECATE_galaxies_hecatecut.parquet')
    assignments = assign_frbs_to_hosts(
        frb_df          = frbs,
        galaxy_catalog  = galaxies,
        localization    = (25., 2., 12.),
        offset_function = 'exponential',
        scale           = 0.5,
        trim_catalog    = 60 * units.arcmin,
        seed            = 42,
    )

    # Step 3 — run PATH
    path_catalog = load_galaxy_catalog('catalog_dudxmmlss_hecate_DECaL.parquet')
    prior_dict   = {
        'P_O_method': 'inverse',
        'PU':         0.15,
        'theta_PDF':  'exp',
        'scale':      0.5,
        'theta_max':  6.0,
    }
    final_sims = run_path.full(
        frbs       = assignments,
        catalog    = path_catalog,
        prior_dict = prior_dict,
        multi      = True,
        ncpu       = 6,
    )

    # Step 4 — build and save digest
    digest = sim_utils.build_digest(
        raw_sim_results  = final_sims,
        frbs             = frbs,
        hosts            = assignments,
        combined_catalog = galaxies,
        output_fn        = 'digest_chime.parquet',
    )
