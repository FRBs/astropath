****************************************************
Notebook: Simulating FRB Populations from a Radio Survey
****************************************************

This notebook demonstrates how to use ``generate_frbs()`` to simulate
realistic FRB populations for different radio surveys, and how to
visualise the resulting DM, redshift, and host-galaxy magnitude
distributions.  For a full description of the simulation procedure see
Section 3.1 of Andersen+26.

.. note::

   The notebook is located at ``docs/nb/Simulate_Generate_FRBs.ipynb``
   in the astropath repository.

What the notebook covers
========================

**1. Listing available surveys**

Prints the full set of supported survey names from ``SURVEY_GRIDS``,
each of which has its own P(z, DM\ :sub:`EG`) grid provided by the
``frb`` package.

**2. Generating FRBs directly from the P(z, DM) grid**

The simplest use-case: DM\ :sub:`EG` and *z* are drawn jointly from the
survey P(z, DM\ :sub:`EG`) grid, a host absolute magnitude
*M*\ :sub:`r` is sampled from the FRB host luminosity function
(``Lz_host_data.csv``), and an apparent magnitude *m*\ :sub:`r` is
computed from the distance modulus at the sampled redshift.

**3. Comparing CHIME, DSA, and CRAFT ICS (10,000 FRBs per survey)**

Side-by-side figures comparing:

* Redshift distributions
* DM\ :sub:`EG` distributions
* DM\ :sub:`EG` vs *z* scatter
* Host apparent magnitude distributions, with Pan-STARRS (~22 mag) and
  DECaLs (~24 mag) survey depth limits overlaid

**4. Generating FRBs from an observed DM catalog**

When simulating a specific observational subsample (e.g. CHIME Catalog 1
bursts with S/N > 12 that trigger baseband saving), a Gaussian KDE is
fitted to the observed DM\ :sub:`EG` values and used to draw synthetic
DMs; redshifts are then drawn from P(z | DM\ :sub:`EG`).  A
side-by-side comparison with the full-grid default is shown, with the
KDE overlaid.

**5. Custom cosmology**

Apparent magnitudes computed with Planck18 vs WMAP9 are compared to
quantify the effect of cosmology choice on *m*\ :sub:`r`.

**6. Saving the output**

Example commands for writing the output DataFrame to CSV and parquet.

Prerequisites
=============

* ``astropath`` installed with the ``frb`` package dependency
* No catalog downloads required for this notebook

Key functions used
==================

* ``generate_frbs(n_frbs, survey, dm_catalog=None, cosmo=None, seed=None)``
* ``SURVEY_GRIDS`` — dictionary mapping survey names to P(z, DM\ :sub:`EG`)
  grid filenames

See :doc:`simulations` for the full API reference and output column
descriptions.
