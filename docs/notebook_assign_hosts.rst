*****************************************************
Notebook: Assigning Simulated FRBs to Host Galaxies
*****************************************************

This notebook shows how to place a simulated FRB population on-sky by
assigning each FRB to a host galaxy from a wide-field optical catalog,
applying a chosen intrinsic galactocentric offset distribution, and
adding a localization error ellipse to produce the observed FRB
position.  For a full description see Section 3.2 of Andersen+26.

.. note::

   The notebook is located at ``docs/nb/Simulate_Assign_Hosts.ipynb``
   in the astropath repository.

What the notebook covers
========================

**Step 1 — Generate an FRB population**

10,000 CHIME FRBs are generated using the DM\ :sub:`EG` distribution of
the non-repeating CHIME/FRB Catalog 2 bursts with S/N > 12 (with a 30
pc cm\ :sup:`−3` Milky Way halo contribution subtracted) as the DM
prior, matching the subsample used in Andersen+26.  See
:doc:`notebook_generate_frbs` for details of this step.

**Step 2 — Load the possible host galaxy catalog**

The combined SGA-2020 + Pan-STARRS/DECaLs + HSC-SSP Wide catalog is
loaded with ``load_galaxy_catalog()``.  Two versions are available from
the `Google Drive link <https://drive.google.com/drive/folders/1PKqh8tnDLbtqIuGeoPFEEh60ovs8Zjw8?usp=drive_link>`_
(~150 MB each):

* ``combined_HSC_PS1_SGA_hosts_clean_ps1_iter10.parquet`` — use when
  running PATH with the Pan-STARRS catalog
* ``combined_HSC_DECaLs_SGA_hosts_clean_iter11.parquet`` — use when
  running PATH with the DECaLs catalog (used in the notebook)

The catalog is loaded from ``$FRB_APATH/catalogs/``.  If the real
catalog is unavailable a mock catalog is generated automatically for
illustration purposes.

The notebook also visualises the catalog sky distribution, magnitude
distribution, and angular-size distribution.

**Step 3 — Assign FRBs to host galaxies**

``assign_frbs_to_hosts()`` is called with a CHIME-KKO localization
ellipse (*a* = 25″, *b* = 2″, PA = 12°) and an exponential intrinsic
offset distribution (``scale = 0.5``).  The output DataFrame is
described at the top of the notebook and in :doc:`simulations`.

**Step 4 — Visualising distributions**

Four sets of figures explore how the output depends on the input
parameters:

* *Assigned host properties* — magnitude and angular-size distributions
  of the assigned hosts, compared against the FRB *m*\ :sub:`r`
  distribution.

* *Intrinsic offset distributions* — 2-D (RA vs Dec) and 1-D (radial)
  host-normalised offset scatter for all three
  ``offset_function`` choices (``'exponential'`` with ``scale = 1``,
  and ``'uniform_1d'`` and ``'uniform_2d'`` with ``scale = 6``).

* *Effect of scale* — exponential offsets with scale = 2.0, 1.0, and
  0.5 overlaid, showing how the concentration around galaxy centres
  changes.

* *Effect of localization size* — observed (localization-convolved)
  host-normalised offsets for three localization sizes: CHIME
  full-array (0.1″ × 0.05″), 1″ circular, and CHIME-KKO (25″ × 2″).

**Step 5 — Saving results**

Example commands for writing the assignment DataFrame to CSV and
parquet.

Prerequisites
=============

* ``astropath`` installed
* ``$FRB_APATH`` pointing to the directory whose ``catalogs/``
  subdirectory contains the host-galaxy catalog (or accept the mock
  catalog fallback)
* A completed ``generate_frbs()`` run, or generate FRBs within the
  notebook

Key functions used
==================

* ``generate_frbs(n_frbs, survey, dm_catalog=None, seed=None)``
* ``load_galaxy_catalog(catalog_fn=...)``
* ``assign_frbs_to_hosts(frb_df, galaxy_catalog, localization, offset_function, scale, seed)``

See :doc:`simulations` for the full API reference and output column
descriptions.
