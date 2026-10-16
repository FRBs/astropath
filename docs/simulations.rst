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
survey-specific P(z, DM\ :sub:`EG`) grid.  For each FRB it also draws a
host-galaxy absolute magnitude *M*\ :sub:`r` from the luminosity
distribution of published FRB hosts, then converts to apparent
magnitude *m*\ :sub:`r` using the sampled redshift and a cosmological
distance modulus.

The host absolute magnitudes are drawn from a Gaussian KDE
(bandwidth 0.3, evaluated over −26 < *M*\ :sub:`r` < −13) of the
published host sample in ``astropath/data/simulations/published_frb_hosts_100726.csv``
(the Dong+26 sample), whose apparent *r*-band magnitudes are converted to
absolute magnitudes with Planck18 luminosity distances.  See Section 3.1
of Andersen+26.

Supported surveys
-----------------

For ``CHIME`` the grid is the James+23 best-fit P(z, DM\ :sub:`EG`)
grid shipped with astropath
(``astropath/data/simulations/CHIME_progenitors_james23_case_bestfit.npz``).
For every other survey the grid is taken from the ``frb`` package
(``frb.dm.prob_dmz``).

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

``generate_frbs(n_frbs, survey, dm_catalog=None, cosmo=None, seed=None)``

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
       drawn from a Gaussian KDE fitted to this array (in
       log\ :sub:`10`\ (DM + 60) space) and redshifts are then drawn from
       P(z | DM\ :sub:`EG`).  When ``None``, DM and *z* are drawn jointly
       from the P(DM, z) grid *(np.ndarray, optional)*
   * - *cosmo*
     - Planck18
     - Astropy cosmology for the distance modulus used to compute
       *m*\ :sub:`r` *(astropy.cosmology, optional)*
   * - *seed*
     - ``None``
     - Global random seed for reproducibility *(int, optional)*

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

Constrain the DM distribution to an observed sample (e.g. the
non-repeating CHIME/FRB Catalog 2 bursts with S/N > 12, with a 30
pc cm\ :sup:`−3` Milky Way halo contribution subtracted)::

    from importlib.resources import files
    import pandas as pd
    import numpy as np

    fn = files('astropath.data') / 'frb_surveys' / 'chimefrbcat2.csv'
    df_dr2 = pd.read_csv(fn)
    df_dr2 = df_dr2[(df_dr2['bonsai_snr'] > 12.) & df_dr2['repeater_name'].isna()]
    dm_halo = 30.  # pc/cc
    dms_eg = np.nanmean([df_dr2['dm_exc_ne2001'].values,
                         df_dr2['dm_exc_ymw16'].values], axis=0) - dm_halo
    observed_dms = dms_eg[dms_eg > 0.]

    frbs = generate_frbs(1000, 'CHIME', dm_catalog=observed_dms, seed=42)

The same DM selection is available as
``astropath.simulations.generate_frbs.load_chime_cat2_DMeg(cut_snr=12., dm_halo=30.)``.


.. _step2:

Step 2 — Load a Galaxy Catalog and Assign FRBs to Hosts
=========================================================

Loading the host galaxy catalog
--------------------------------

``load_galaxy_catalog()`` reads a parquet catalog from the ``catalogs``
subdirectory of the directory given by the ``FRB_APATH`` environment
variable, i.e. ``$FRB_APATH/catalogs/<catalog_fn>``.  ``catalog_fn``
defaults to ``combined_HSC_DECaLs_SGA_hosts_clean_iter11.parquet``.  It
raises a ``ValueError`` if ``FRB_APATH`` is not set or the file does not
exist::

    import os
    from astropath.simulations import load_galaxy_catalog

    os.environ['FRB_APATH'] = '/path/to/simulation/data/'
    galaxies = load_galaxy_catalog(
        catalog_fn='combined_HSC_DECaLs_SGA_hosts_clean_iter11.parquet'
    )

Two combined host catalogs are available from the
`Google Drive link <https://drive.google.com/drive/folders/1PKqh8tnDLbtqIuGeoPFEEh60ovs8Zjw8?usp=drive_link>`_
(~150 MB each).  Choose the one that matches the PATH catalog you will
use in Step 3:

* ``combined_HSC_PS1_SGA_hosts_clean_ps1_iter10.parquet`` — pair with
  the Pan-STARRS PATH catalog
* ``combined_HSC_DECaLs_SGA_hosts_clean_iter11.parquet`` — pair with
  the DECaLs PATH catalog

Both catalogs combine SGA-2020 large galaxies, Pan-STARRS or DECaLs
galaxies, and HSC-SSP Wide galaxies.  See Section 3.2 of Andersen+26
for how they are constructed.

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
apparent magnitude, places the FRB within the host according to a
chosen intrinsic offset distribution, and offsets the position by a
localization error to produce the *observed* FRB coordinates.

Magnitude matching is a single pass without replacement: each FRB takes
the unused galaxy closest to its *m*\ :sub:`r`, with the most
oversubscribed magnitudes (highest ratio of FRBs to galaxies per 0.25
mag) served first.  Where the catalog cannot supply an unused galaxy
within 0.05 mag of the ideal match, the FRB re-uses the least-used
galaxy within that tolerance.

The localization error is drawn from an elliptical Gaussian with 1σ
widths *a* along the position angle PA and *b* along PA + 90°.

API
~~~

``assign_frbs_to_hosts(frb_df, galaxy_catalog, localization, mag_range=None, offset_function='exponential', scale=0.5, dm_neighbors=25, dm_match_scale='log', seed=None, debug=False)``

.. list-table::
   :header-rows: 1
   :widths: 22 18 60

   * - Parameter
     - Default
     - Description
   * - *frb_df*
     - —
     - FRB DataFrame from ``generate_frbs()``; must contain column
       ``m_r`` (and ``DMeg`` for joint DM–localization sampling)
       *(pd.DataFrame)*
   * - *galaxy_catalog*
     - —
     - Galaxy catalog DataFrame; see column requirements above
       *(pd.DataFrame)*
   * - *localization*
     - —
     - Error ellipse(s), in one of three forms: a single ``(a, b, PA)``
       triple applied to every FRB, where *a* and *b* are the semi-major
       and semi-minor axes in arcsec and *PA* is the position angle in
       degrees East of North; a population of ellipses (an (N, 3) array or
       a DataFrame with columns ``a``, ``b``, ``PA``), from which one row
       is drawn at random per FRB; or a joint (DM, a, b, PA) population
       (an (N, 4) array, or a DataFrame that also has a ``DM`` column), in
       which each FRB is given the localization of a real FRB of similar
       DM.  ``localizations_from_errors()`` builds either population from
       an observed catalog *(tuple, array, or pd.DataFrame)*
   * - *mag_range*
     - ``None``
     - ``(min, max)`` apparent magnitude range; FRBs outside this range
       are excluded.  ``None`` keeps all FRBs *(tuple, optional)*
   * - *offset_function*
     - ``'exponential'``
     - Intrinsic galactocentric offset distribution.
       One of ``'exponential'``, ``'exponential_incorrect'``,
       ``'uniform_1d'``, ``'uniform_2d'``; see below *(str)*
   * - *scale*
     - ``0.5``
     - Scale parameter for the chosen offset distribution, in units of
       the host half-light radius *(float)*
   * - *dm_neighbors*
     - ``25``
     - Joint DM–localization sampling only: number of real FRBs, nearest
       in DM, from which each localization is drawn *(int)*
   * - *dm_match_scale*
     - ``'log'``
     - Joint DM–localization sampling only: ``'log'`` or ``'linear'``,
       the space in which DM neighbours are found *(str)*
   * - *seed*
     - ``None``
     - Random seed *(int, optional)*
   * - *debug*
     - ``False``
     - Print matching and localization diagnostics *(bool)*

Intrinsic offset distributions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - ``offset_function``
     - Distribution
   * - ``'exponential'``
     - :math:`p(r) \propto r\,\exp(-r/s)` — Gamma(2, scale) in units of
       the host half-light radius, truncated at 6 half-light radii.
       Matches the PATH exponential prior.
   * - ``'exponential_incorrect'``
     - :math:`p(r) \propto \exp(-r/s)` — exponential in *r* in units of
       the host half-light radius, truncated at 6 half-light radii.  Does
       not match the PATH exponential prior.
   * - ``'uniform_1d'``
     - Uniform in *r*, truncated at *scale* × half-light radius (at most
       10 half-light radii).
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
     - Localization semi-major axis used for this FRB (arcsec)
   * - ``b``
     - float
     - Localization semi-minor axis used for this FRB (arcsec)
   * - ``PA``
     - float
     - Localization position angle used for this FRB (degrees)

Usage example::

    from astropath.simulations import assign_frbs_to_hosts

    assignments = assign_frbs_to_hosts(
        frb_df          = frbs,
        galaxy_catalog  = galaxies,
        localization    = (25., 2., 12.),   # a, b, PA  [arcsec, arcsec, deg]
        offset_function = 'exponential',
        scale           = 0.5,
        seed            = 42,
    )


.. _step3:

Step 3 — Run PATH
==================

``run_path.full()`` runs PATH on every simulated FRB scenario.  The
computation is embarrassingly parallel and is distributed across CPUs
with Python's ``multiprocessing`` module.  As a reference point: 10,000
FRBs with a 50 × 100 mas localization took ~2 minutes on an Apple M3
MacBook Pro (16 GB RAM, ``ncpu = 6``), and ~6 minutes for a 25″ × 25″
localization.  Memory usage also increases rapidly with localization size
and number of FRBs.  For example, 100,000 FRBs processing on 10 CPU cores
on CANFAR can take >100 GB of RAM.

Each FRB is run in PATH's ``'local'`` likelihood mode with a relative
step size of 0.05, against the PATH-catalog galaxies within a search
radius of max(10 *a*, 3′) of its observed position.

Before calling this function, load the pre-queried PATH galaxy catalog.
Two options are available from the
`Google Drive link <https://drive.google.com/drive/folders/1PKqh8tnDLbtqIuGeoPFEEh60ovs8Zjw8?usp=drive_link>`_:

* ``catalog_hscwide_SGA_PS1_clean_ps1_iter10.parquet`` (Pan-STARRS)
* ``catalog_hscwide_SGA_DECaL_clean_iter11.parquet`` (DECaLs)

Choose the catalog that matches the host-galaxy catalog used in Step 2,
save it to ``$FRB_APATH/catalogs/``, and load it::

    path_catalog = load_galaxy_catalog(
        catalog_fn='catalog_hscwide_SGA_DECaL_clean_iter11.parquet'
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

``run_path.full(frbs, catalog, prior_dict, multi=True, ncpu=4, frb_chunk=5000, debug=False)``

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
   * - *frb_chunk*
     - ``5000``
     - Number of FRBs cross-matched against the catalog at a time when
       building each FRB's candidate list *(int)*
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
     - ``'inverse'``, ``'inverse_diff'``, ``'identical'``,
       ``'updated_driver'``, ``'updated_driver_diff'``, ``'user'``
     - Magnitude prior :math:`P(O_i)`.  ``'inverse'`` weights by the
       inverse of the cumulative galaxy surface density brighter than the
       candidate; ``'inverse_diff'`` by the inverse of the differential
       galaxy counts at the candidate's magnitude; ``'identical'`` gives equal weight to all candidates;
       ``'updated_driver'`` and ``'updated_driver_diff'`` use a fit to the
       Driver+16 r-band counts (cumulative and differential,
       respectively); ``'user'`` calls the user-defined function
       ``astropath.priors.USR_raw_prior_Oi``.
   * - ``'PU'``
     - float in [0, 1)
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

Each row is one PATH candidate for one simulated FRB.  The candidates
of each FRB are sorted by descending :math:`P(O_i | x)`, so the first
row for each ``iFRB`` is the best candidate.  Key columns:

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
   * - ``ID``
     - int
     - ID of the candidate galaxy in the PATH catalog
   * - ``sep``
     - float
     - Separation between the candidate and the observed FRB position
       (arcsec)
   * - ``P_O``
     - float
     - PATH prior :math:`P(O_i)`
   * - ``P_Ox``
     - float
     - PATH posterior :math:`P(O_i | x)`
   * - ``P_Ux``
     - float
     - PATH unseen-host posterior :math:`P(U | x)`
   * - ``iFRB``
     - int
     - Row position (0-based) of the FRB in the ``frbs`` DataFrame

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
into a single, analysis-ready DataFrame with one row per simulated FRB
(its best PATH candidate), and optionally writes it to a parquet file.
FRBs with no PATH candidates are excluded.

API
---

``utils.build_digest(raw_sim_results, frbs, hosts, combined_catalog, output_fn=None)``

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

An association is "correct" when the best candidate is the true host:
by ID where the true host's ID appears among the candidates PATH
searched (for any FRB in the run), and otherwise when the two are separated by less than
max(*φ*\ :sub:`host`, *φ*\ :sub:`cand`, 0.5″), the larger of the two
half-light radii with a 0.5″ floor.  A host is "unseen" when its ID is
not among the searched candidates and the nearest searched candidate
lies farther than max(*φ*\ :sub:`host`, *φ*\ :sub:`nearest`, 0.5″).
An unseen host is never a correct association.

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
   * - ``sep_host_loc_norm``
     - ``sep_host_loc_arcsec`` normalised by ``ang_size_host``
   * - ``sep_best_loc_arcsec``
     - Separation between the best candidate center and the localization (arcsec)
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
   * - ``ang_size_cand``
     - Angular size of the best candidate (arcsec)
   * - ``mag_cand``
     - Apparent magnitude of the best candidate
   * - ``cand_ID``
     - ID of the best candidate in the PATH catalog
   * - ``P_O``
     - PATH prior :math:`P(O_i)` for the best candidate
   * - ``P_Ox``
     - PATH posterior :math:`P(O_i | x)` for the best candidate
   * - ``P_Ux``
     - PATH unseen-host posterior :math:`P(U | x)`
   * - ``unseen``
     - ``True`` if the true host is not detected in the catalog PATH
       searched (see above)
   * - ``correct_association``
     - ``True`` if the best candidate is the true host (see above)
   * - ``host_catalog``
     - Catalog the true host was drawn from: ``'SGA'``, ``'DECaL_ID'``,
       ``'Pan-STARRS_ID'``, ``'HSC_ID'``, or ``'HECATE_ID'``
       (``'No Match'`` if the host ID is found in none of them, and
       ``'Random'`` for FRBs placed at random positions)

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

``astropath/simulations/run_path_simulation_cli_final.py`` wraps the
entire four-step pipeline into a single command.  It enforces valid
choices and types for every parameter and generates output filenames
automatically from the simulation settings.
``run_path_simulation_cli_final_ps1.py`` is identical except that its
default catalogs are the Pan-STARRS ones.

As for ``load_galaxy_catalog()``, ``FRB_APATH`` must be set in your
environment to the base simulation directory, and every catalog or DM
table named on the command line is read from
``$FRB_APATH/catalogs/<filename>``.  The script exits with an error if
``FRB_APATH`` is not set. If the host or PATH catalog file is not found,
the script prints a warning and falls back to a small mock catalog,
which is only useful for testing. A missing DM or localization/DM table
raises an error.

Usage
-----

.. code-block:: bash

   export FRB_APATH=/path/to/simulation/data   # contains catalogs/

   python run_path_simulation_cli_final.py \
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

The DMs are taken from, in order of precedence:

* ``--locdm-catalog`` — DMs *and* localizations from a table of real
  FRBs (the localization of each simulated FRB is drawn from real FRBs of
  similar DM; ``--loc-a/b/pa`` are ignored);
* ``--dm-catalog`` — DMs only; the localization is ``--loc-a/b/pa``;
* otherwise, for ``--survey CHIME``, the non-repeating CHIME/FRB Catalog 2
  bursts with S/N > 12 (see :ref:`step1`);
* otherwise, the survey's full P(DM, z) grid.

``--dm-catalog`` and ``--locdm-catalog`` cannot be combined.

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
   * - ``--locdm-catalog``
     - ``None``
     - Parquet or CSV table of real FRBs with columns ``DM_ex``
       (extragalactic DM), ``a_err``, ``b_err`` (arcsec), and ``pa``
       (deg).  DMs are drawn from its ``DM_ex`` values and localizations
       are sampled jointly with DM from its ellipses, overriding
       ``--loc-a/b/pa``
   * - ``--dm-catalog``
     - ``None``
     - Parquet or CSV table with an extragalactic-DM column ``DM_ex``.
       DMs are drawn from these values; the localization is
       ``--loc-a/b/pa``
   * - ``--offset-dist``
     - ``exponential``
     - Intrinsic offset distribution (``exponential``,
       ``exponential_incorrect``, ``uniform_1d``, ``uniform_2d``)
   * - ``--offset-scale``
     - 0.5
     - Scale of the intrinsic offset distribution
   * - ``--mag-prior``
     - ``inverse``
     - PATH magnitude prior (``inverse``, ``inverse_diff``,
       ``identical``, ``user``, ``updated_driver``,
       ``updated_driver_diff``).  ``user`` divides a KDE of the
       simulated FRB host *m*\ :sub:`r` distribution by a fit to the
       Driver+16 differential galaxy counts
   * - ``--unseen-prior``
     - 0.15
     - PATH unseen-host prior P(U); must be in [0, 1]
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
     - Filename of the wide-field host catalog
   * - ``--path-catalog``
     - (see note)
     - Filename of the pre-queried PATH galaxy catalog
   * - ``--output-dir``
     - ``./``
     - Output directory
   * - ``--full-output``
     - auto
     - Parquet path for the raw PATH results
   * - ``--digest-output``
     - auto
     - Parquet path for the digest
   * - ``--tag``
     - ``sim``
     - Label placed at the start of the auto-generated filenames
   * - ``--save-full-output``
     - ``True``
     - Whether to write the raw PATH results (can be large); use
       ``--save-full-output False`` to skip
   * - ``--seed``
     - ``None``
     - Global random seed
   * - ``--no-multi``
     - off
     - Disable multiprocessing (serial execution)
   * - ``--debug``
     - off
     - Process only 100 FRBs; useful for testing
   * - ``--random-locs``
     - ``False``
     - Place FRBs at random positions in the XMM-LSS HSC-SSP field
       instead of in host galaxies; ``--random-locs`` or
       ``--random-locs True`` turns it on

Default catalog filenames (looked up in ``$FRB_APATH/catalogs/``):

* ``run_path_simulation_cli_final.py``: host catalog
  ``combined_HSC_DECaLs_SGA_hosts_clean_iter11.parquet``, PATH catalog
  ``catalog_hscwide_SGA_DECaL_clean_iter11.parquet``
* ``run_path_simulation_cli_final_ps1.py``: host catalog
  ``combined_HSC_PS1_SGA_hosts_clean_ps1_iter10.parquet``, PATH catalog
  ``catalog_hscwide_SGA_PS1_clean_ps1_iter10.parquet``

The suite scripts (:ref:`full_suite`) pass ``--tag
<table1_label>_<unseen_prior>``, so their auto-generated output filenames
follow the pattern::

    digest_<table1_label>_<unseen_prior>_<loc>_<intrinsic_offset_distribution>_<offset_prior>_<magnitude_prior>_<path_catalog>_<NFRBs>_<random_seed>.parquet
    full_path_results_<table1_label>_<unseen_prior>_<loc>_<intrinsic_offset_distribution>_<offset_prior>_<magnitude_prior>_<path_catalog>_<NFRBs>_<random_seed>.parquet

where:

* ``<table1_label>`` is the simulation label from Table 1 of Andersen+26
  (e.g. ``default-CK``) and ``<unseen_prior>`` is ``pu<P(U)>`` (e.g.
  ``pu0.10``).  Together they form ``--tag``; a run with a different
  ``--tag`` starts its filenames with that tag instead.  P(U) appears in
  the filename only through ``--tag``.
* ``<loc>`` is ``a<a>b<b>pa<PA>`` built from ``--loc-a/b/pa`` with two
  decimals each (also when ``--locdm-catalog`` is used).
* ``<intrinsic_offset_distribution>``, ``<offset_prior>`` and
  ``<magnitude_prior>`` are ``--offset-dist``, ``--offset-prior`` and
  ``--mag-prior``.
* ``<path_catalog>`` is the ``--path-catalog`` filename without
  ``.parquet``, ``<NFRBs>`` is ``--n-frbs``, and ``<random_seed>`` is
  ``seed<seed>`` (e.g. ``seed57822``).

Run ``python run_path_simulation_cli_final.py --help`` for the complete
argument list.


.. _andersen26:

Reproducing Andersen+26 Figures
================================

The notebook ``docs/nb/Reproducing_Andersen+26_Figures.ipynb`` contains
code to reproduce Figures 3–16 of Andersen+26.  To run it, download the
pre-computed simulation digests from the
`Google Drive link <https://drive.google.com/drive/folders/1Rv6koYfKJ7yQ5xYV666gIaMQ2wX0A6Tg?usp=drive_link>`_
and save them to ``$FRB_APATH/sim_paper/``.  Digest files are named::

    digest_<table1_label>_<unseen_prior>_<loc>_<intrinsic_offset_distribution>_<offset_prior>_<magnitude_prior>_<path_catalog>_<NFRBs>_<random_seed>.parquet

For example,
``digest_default-CK_pu0.10_a25.00b2.00pa12.00_exponential_exp_inverse_catalog_hscwide_SGA_DECaL_clean_iter11_100000_seed57822.parquet``
is the ``default-CK`` simulation from Table 1 of Andersen+26: the
CHIME-KKO localization, P(U) = 0.10, an exponential intrinsic offset
distribution, an exponential PATH offset prior, and an inverse magnitude
prior, run with the DECaLs PATH catalog.

Some figures need further files:

* Figures 4 and 5 read the DECaLs host catalog, Figure 5 also the DECaLs
  PATH catalog, and Figure 7 the DECaLs and Pan-STARRS PATH catalogs, all
  from ``$FRB_APATH/catalogs/``; see :ref:`step2` and :ref:`step3` for
  download instructions.  The PATH catalogs are large: 1.1 GB for
  Pan-STARRS and 5.1 GB for DECaLs.
* Figures 10, 13, and 16 read the raw PATH results
  (``full_path_results_<...>.parquet``, named like the digests) of the
  ``default-CK`` and ``mismatch-CKGH`` simulations from
  ``$FRB_APATH/sim_paper/``.
* Figure 7 reads the generated DSA and CRAFT FRB populations
  ``generated_frbs_dsa_redo.parquet`` and
  ``generated_frbs_craft_redo.parquet`` from ``$FRB_APATH/sim_paper/``.

Figures are saved to ``$FRB_APATH/paper_figures/``, which must exist.


Complete Example
================

The following reproduces a typical CHIME simulation end-to-end::

    import os
    import numpy as np
    import pandas as pd
    from importlib.resources import files
    from astropath.simulations import generate_frbs, assign_frbs_to_hosts
    from astropath.simulations import load_galaxy_catalog, run_path
    from astropath.simulations import utils as sim_utils

    os.environ['FRB_APATH'] = '/path/to/simulation/data/'

    # Step 1 — generate FRBs
    fn = files('astropath.data') / 'frb_surveys' / 'chimefrbcat2.csv'
    df_dr2 = pd.read_csv(fn)
    df_dr2 = df_dr2[(df_dr2['bonsai_snr'] > 12.) & df_dr2['repeater_name'].isna()]
    dms_eg = np.nanmean([df_dr2['dm_exc_ne2001'].values,
                         df_dr2['dm_exc_ymw16'].values], axis=0) - 30.
    observed_dms = dms_eg[dms_eg > 0.]
    frbs = generate_frbs(1000, 'CHIME', dm_catalog=observed_dms, seed=42)

    # Step 2 — assign to hosts
    galaxies    = load_galaxy_catalog(
                      'combined_HSC_DECaLs_SGA_hosts_clean_iter11.parquet')
    assignments = assign_frbs_to_hosts(
        frb_df          = frbs,
        galaxy_catalog  = galaxies,
        localization    = (25., 2., 12.),
        offset_function = 'exponential',
        scale           = 0.5,
        seed            = 42,
    )

    # Step 3 — run PATH
    path_catalog = load_galaxy_catalog('catalog_hscwide_SGA_DECaL_clean_iter11.parquet')
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


.. _full_suite:

Running the Full Simulation Suite
=================================

The simulations in Table 1 of Andersen+26 are run with two shell
scripts in ``astropath/simulations/``, each of which calls the CLI
(:ref:`cli`) once per simulation:

* ``run_all_simulations_final.sh`` — DECaLs runs, calling
  ``run_path_simulation_cli_final.py`` (seed 57822, P(U) = 0.10)
* ``run_all_simulations_final_ps1.sh`` — Pan-STARRS runs, calling
  ``run_path_simulation_cli_final_ps1.py`` (seed 57821, P(U) = 0.30)

Both scripts run every simulation for ``CHIME`` with 100,000 FRBs on 12
CPUs, an intrinsic offset scale of 0.5, and an ``exp`` PATH offset
prior with ``scale = 0.5`` and ``theta_max = 6.0``.  The localizations
are:

.. list-table::
   :header-rows: 1
   :widths: 20 30 50

   * - Label
     - (*a*, *b*, PA)
     - Description
   * - ``C``
     - (25″, 25″, 0°)
     - CHIME-only
   * - ``CK``
     - (25″, 2″, 12°)
     - CHIME-KKO
   * - ``CKG``
     - (15″, 0.05″, 5°)
     - CHIME-KKO-GBO
   * - ``CKGH``
     - (0.1″, 0.05″, 10°)
     - CHIME-KKO-GBO-HCRO
   * - ``1arcsec``
     - (1″, 1″, 0°)
     - Generic 1″

The simulations are organised in groups, each a block of ``run_sim``
calls:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Group
     - Simulations
   * - ``default-*``
     - All five localizations, ``inverse`` magnitude prior
   * - ``user-*``
     - All five localizations, ``user`` magnitude prior (the simulated
       host *m*\ :sub:`r` distribution; see :ref:`cli`)
   * - ``updateddriver-*``
     - ``C`` and ``CK`` with the ``updated_driver`` and
       ``updated_driver_diff`` priors (DECaLs script only)
   * - ``mismatch-CKGH``, ``ident-1arcsec``
     - ``CKGH`` with a ``uniform_1d`` intrinsic offset distribution
       (scale 8.0) against the ``exp`` prior; ``1arcsec`` with the
       ``identical`` magnitude prior
   * - ``ident-*``
     - All five localizations, ``identical`` magnitude prior
   * - ``PU-*``
     - All five localizations, ``inverse`` prior, P(U) swept from 0.0 to
       0.9 in steps of 0.1, with 10,000 FRBs per run

Each simulation's label (e.g. ``default-CK_pu0.10``) is passed as
``--tag``, so it heads the output filenames (see :ref:`andersen26`).

To run the suite:

1. In the configuration block at the top of the shell script, set
   ``FRB_APATH`` (exported to the CLI) to your base simulation directory,
   the one containing ``catalogs/``.  Outputs are written to
   ``OUTPUT_DIR = $FRB_APATH/sim_paper``, where the Reproducing notebook
   expects them.
2. Edit the rest of the configuration block if needed (``N_FRBS``, ``NCPU``, ``SEED``, the global ``pu``, and
   ``PU_VALUES`` for the ``PU-*`` sweep).
3. Comment out the groups you do not want to run.  In
   ``run_all_simulations_final.sh`` every group is active; in
   ``run_all_simulations_final_ps1.sh`` only ``default-*`` is active and
   the other groups are commented out.  ``TOTAL`` only sets the progress
   counter printed for each run.
4. Run the script from ``astropath/simulations/``, or set ``CLI`` to the
   path of the CLI script::

       cd astropath/simulations
       bash run_all_simulations_final.sh

       # or, from elsewhere
       CLI=/path/to/run_path_simulation_cli_final.py bash run_all_simulations_final.sh

The script stops at the first failed run (``set -euo pipefail``).  Each
run writes a digest and, by default, the raw PATH results to
``OUTPUT_DIR``.
