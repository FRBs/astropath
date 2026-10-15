*****************************************
Notebook: Running PATH on Simulated FRBs
*****************************************

This notebook walks through the complete end-to-end PATH simulation
pipeline — from generating FRBs to evaluating association performance
— and is the primary reference notebook for reproducing the simulations
in Andersen+26.

.. note::

   The notebook is located at ``docs/nb/Simulate_Run_PATH.ipynb``
   in the astropath repository.

What the notebook covers
========================

**Step 1 — Generate FRB population**

5,000 CHIME FRBs are generated using the DM\ :sub:`EG` distribution of
the non-repeating CHIME/FRB Catalog 2 bursts with S/N > 12 (with a 30
pc cm\ :sup:`−3` Milky Way halo contribution subtracted) as the DM
prior.  See :doc:`notebook_generate_frbs` for details.

**Step 2 — Load the host galaxy catalog**

The combined SGA/DECaLs/HSC host catalog
(``combined_HSC_DECaLs_SGA_hosts_clean_iter11.parquet``) is loaded from
``$FRB_APATH/catalogs/`` via ``load_galaxy_catalog()``, with a
mock-catalog fallback for demonstration.  See
:doc:`notebook_assign_hosts` and Section 3.2 of Andersen+26 for
details.

**Step 3 — Assign FRBs to host galaxies**

FRBs are assigned to hosts with a CHIME-KKO localization
(*a* = 25″, *b* = 2″, PA = 12°) and an exponential intrinsic
offset distribution (``scale = 0.5``).  See
:doc:`notebook_assign_hosts` and Section 3.2 of Andersen+26 for
details.

**Step 4 — Run PATH**

The pre-queried PATH galaxy catalog must be downloaded from the
`Google Drive link <https://drive.google.com/drive/folders/1PKqh8tnDLbtqIuGeoPFEEh60ovs8Zjw8?usp=drive_link>`_
and saved to ``$FRB_APATH/catalogs/``:

* ``catalog_hscwide_SGA_PS1_clean_ps1_iter10.parquet`` (Pan-STARRS)
* ``catalog_hscwide_SGA_DECaL_clean_iter11.parquet`` (DECaLs; used in
  the notebook)

``run_path.full()`` is called with the prior dictionary shown below.
Each simulated FRB is an independent PATH run so the workload is
parallelised with Python's ``multiprocessing`` module.  As a reference:
5,000 FRBs took ~50 minutes on an Apple M3 MacBook Pro (16 GB RAM,
``ncpu = 6``).

.. code-block:: python

   prior_dict = {
       'P_O_method': 'inverse',
       'PU':         0.15,
       'theta_PDF':  'exp',
       'scale':      0.5,
       'theta_max':  6.0,
   }

See Section 3.3 of Andersen+26 for details on prior selection.

**Step 5 — Build the simulation digest**

``sim_utils.build_digest()`` merges all pipeline DataFrames into a
single analysis-ready table, which is then saved as
``simulation_digest.parquet``.  The digest columns are listed at the
top of the notebook and in :doc:`simulations`.

**Step 6 — Analyse simulation results**

The saved digest is re-loaded and used to produce performance plots,
including:

* Correct-association, incorrect-association, and non-association rates
  at a chosen P(O|x) threshold
* PATH posterior P(O|x) vs true-host and best-candidate apparent
  magnitude (split by correct vs incorrect associations)
* Magnitude bias analysis
* Offset bias analysis: the simulated and recovered host-normalised
  offset distributions compared with the PATH offset prior convolved
  with the localization (``sim_utils.stack_convolved_prior()``)
* Posterior calibration: the correct-association fraction vs P(O|x)
  and the unseen fraction vs P(U|x), with bootstrap errors

See ``Reproducing_Andersen+26_Figures.ipynb`` (:ref:`andersen26`) for
more detailed analysis.

Prerequisites
=============

* ``astropath`` installed with the ``frb`` package dependency
* ``$FRB_APATH`` pointing to the directory whose ``catalogs/``
  subdirectory contains both the host-galaxy catalog and the
  pre-queried PATH catalog
* Sufficient RAM and CPU for the PATH step (~16 GB RAM recommended for
  large runs)

Key functions used
==================

* ``generate_frbs(n_frbs, survey, dm_catalog=None, seed=None)``
* ``load_galaxy_catalog(catalog_fn=...)``
* ``assign_frbs_to_hosts(frb_df, galaxy_catalog, localization, offset_function, scale, seed)``
* ``run_path.full(frbs, catalog, prior_dict, multi, ncpu)``
* ``sim_utils.build_digest(raw_sim_results, frbs, hosts, combined_catalog)``
* ``sim_utils.stack_convolved_prior(u, theta_prior, frb_df, phi_col)``

See :doc:`simulations` for the full API reference, prior dictionary
options, and digest column descriptions.
