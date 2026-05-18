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

1,000 CHIME FRBs are generated using the CHIME Catalog 1
DM\ :sub:`EG` distribution (S/N > 12) as the DM prior.  See
:doc:`notebook_generate_frbs` for details.

**Step 2 — Load the host galaxy catalog**

The combined HSC/DECaLs/HECATE catalog is loaded via
``load_galaxy_catalog()``, with a mock-catalog fallback for
demonstration.  See :doc:`notebook_assign_hosts` and Section 3.2 of
Andersen+26 for details.

**Step 3 — Assign FRBs to host galaxies**

FRBs are assigned to hosts with a CHIME-KKO localization
(*a* = 25″, *b* = 2″, PA = 12°) and an exponential intrinsic
offset distribution (``scale = 0.5``).  See
:doc:`notebook_assign_hosts` and Section 3.2 of Andersen+26 for
details.

**Step 4 — Run PATH**

The pre-queried PATH galaxy catalog must be downloaded from the
`Google Drive link <https://drive.google.com/drive/folders/1PKqh8tnDLbtqIuGeoPFEEh60ovs8Zjw8?usp=drive_link>`_
and saved to ``$FRB_APATH``:

* ``catalog_dudxmmlss_hecate_Pan-STARRS.parquet`` (900 MB)
* ``catalog_dudxmmlss_hecate_DECaL.parquet`` (6.1 GB)

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
single analysis-ready table and saves it as
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

Prerequisites
=============

* ``astropath`` installed with the ``frb`` package dependency
* ``$FRB_APATH`` pointing to the directory containing both the
  host-galaxy catalog and the pre-queried PATH catalog
* Sufficient RAM and CPU for the PATH step (~16 GB RAM recommended for
  large runs)

Key functions used
==================

* ``generate_frbs(n_frbs, survey, dm_catalog=None, seed=None)``
* ``load_galaxy_catalog(catalog_fn=...)``
* ``assign_frbs_to_hosts(frb_df, galaxy_catalog, localization, offset_function, scale, trim_catalog, seed)``
* ``run_path.full(frbs, catalog, prior_dict, multi, ncpu)``
* ``sim_utils.build_digest(raw_sim_results, frbs, hosts, combined_catalog, output_fn)``

See :doc:`simulations` for the full API reference, prior dictionary
options, and digest column descriptions.
