***********
Simulations
***********

This document describes the simulation tools available in *astropath*
for generating synthetic FRB populations.

Overview
========

The ``astropath.simulations`` module provides tools for generating
simulated Fast Radio Burst (FRB) populations with realistic properties.
This is useful for:

* Testing PATH analysis pipelines
* Estimating detection efficiencies
* Studying selection effects
* Planning observations

The primary function is ``generate_frbs()`` which produces FRBs with
the following properties:

* **DM**: Extragalactic dispersion measure (pc/cm\ :sup:`3`)
* **z**: Redshift
* **M_r**: Host galaxy absolute r-band magnitude
* **m_r**: Host galaxy apparent r-band magnitude

Supported Surveys
=================

The module supports generation of FRBs for the following surveys:

* **CHIME**: Canadian Hydrogen Intensity Mapping Experiment
* **DSA**: Deep Synoptic Array (DSA-110)
* **ASKAP**: Australian SKA Pathfinder (CRAFT surveys)
* **CRAFT**: Alias for ASKAP
* **CRAFT_ICS_1300**: CRAFT Incoherent Sum at 1300 MHz
* **CRAFT_ICS_892**: CRAFT Incoherent Sum at 892 MHz
* **CRAFT_ICS_1632**: CRAFT Incoherent Sum at 1632 MHz
* **Parkes**: Parkes Multibeam
* **FAST**: Five-hundred-meter Aperture Spherical Telescope

Each survey has a unique P(z,DM) grid that captures the selection
effects and sensitivity of that instrument.

Basic Usage
===========

The simplest usage generates FRBs by sampling directly from the
survey-specific P(z,DM) grid::

    from astropath.simulations import generate_frbs

    # Generate 1000 CHIME FRBs
    df = generate_frbs(1000, 'CHIME', seed=42)

    # View the results
    print(df.head())
    #           DM         z       M_r        m_r
    # 0   234.0000  0.150000 -19.85432  18.234521
    # 1   567.0000  0.420000 -21.23456  21.876543
    # ...

The output is a pandas DataFrame with columns for DM, z, M_r, and m_r.

Generating for Different Surveys
================================

Simply change the survey name to generate FRBs with different
selection functions::

    # DSA-110 FRBs
    df_dsa = generate_frbs(1000, 'DSA', seed=42)

    # ASKAP/CRAFT FRBs
    df_askap = generate_frbs(1000, 'ASKAP', seed=42)

    # Compare redshift distributions
    print(f"CHIME median z: {df['z'].median():.3f}")
    print(f"DSA median z: {df_dsa['z'].median():.3f}")
    print(f"ASKAP median z: {df_askap['z'].median():.3f}")

Using an Observed DM Catalog
============================

If you have observed DM values from a specific sample, you can use
them to constrain the DM distribution via kernel density estimation::

    import numpy as np

    # Example: DMs from a specific observing campaign
    observed_dms = np.array([362.4, 589.0, 364.5, 339.5, 322.2, 594.6])

    # Generate FRBs with DMs sampled from a KDE of the observed values
    df = generate_frbs(1000, 'DSA', dm_catalog=observed_dms, seed=42)

This fits a Gaussian KDE to the provided DM values and samples new
DMs from that distribution, then draws redshifts from the P(z|DM)
grid.

Reproducibility
===============

Use the ``seed`` parameter for reproducible results::

    df1 = generate_frbs(100, 'CHIME', seed=42)
    df2 = generate_frbs(100, 'CHIME', seed=42)

    # These will be identical
    assert (df1['DM'] == df2['DM']).all()

Custom Cosmology
================

By default, simulations use Planck18 cosmology. You can specify
a different cosmology::

    from astropy.cosmology import WMAP9

    df = generate_frbs(1000, 'CHIME', cosmo=WMAP9, seed=42)

Algorithm Details
=================

The ``generate_frbs()`` function follows these steps:

1. **Load P(z,DM) Grid**: Survey-specific grids from the ``frb`` package
   capture the probability of detecting an FRB at redshift z with
   dispersion measure DM, accounting for telescope sensitivity and
   selection effects.

2. **Sample DM Values**: Either directly from the P(DM) marginal
   distribution of the grid, or from a KDE fit to user-provided
   catalog values.

3. **Sample Redshifts**: For each DM, sample z from P(z|DM) using
   inverse transform sampling on the cumulative distribution.

4. **Sample Host M_r**: Draw absolute magnitudes from the observed
   FRB host galaxy luminosity function (loaded from
   ``frb.galaxies.hosts.load_Mr_pdf()``).

5. **Compute Apparent Magnitudes**: Calculate m_r = M_r + distance_modulus(z)
   using the specified cosmology.

Dependencies
============

The simulations module requires the ``frb`` package to be installed::

    pip install frb

This provides access to:

* P(z,DM) grids for various surveys (``frb.dm.prob_dmz``)
* Host galaxy magnitude distributions (``frb.galaxies.hosts``)

API Reference
=============

generate_frbs
-------------

.. py:function:: generate_frbs(n_frbs, survey, dm_catalog=None, cosmo=None, seed=None)

   Generate a population of simulated FRBs.

   :param n_frbs: Number of FRBs to generate
   :type n_frbs: int
   :param survey: Survey name (e.g., 'CHIME', 'DSA', 'ASKAP')
   :type survey: str
   :param dm_catalog: Optional array of observed DMs to sample from
   :type dm_catalog: np.ndarray, optional
   :param cosmo: Cosmology for distance calculations (default: Planck18)
   :type cosmo: astropy.cosmology, optional
   :param seed: Random seed for reproducibility
   :type seed: int, optional
   :returns: DataFrame with columns 'DM', 'z', 'M_r', 'm_r'
   :rtype: pandas.DataFrame

SURVEY_GRIDS
------------

.. py:data:: SURVEY_GRIDS

   Dictionary mapping survey names to their P(z,DM) grid filenames.
   Available surveys: CHIME, DSA, ASKAP, CRAFT, CRAFT_ICS_1300,
   CRAFT_ICS_892, CRAFT_ICS_1632, Parkes, FAST.

Example: Full Workflow
======================

Here is a complete example generating FRBs and examining their
properties::

    import matplotlib.pyplot as plt
    from astropath.simulations import generate_frbs

    # Generate FRBs for three surveys
    surveys = ['CHIME', 'DSA', 'ASKAP']
    dfs = {s: generate_frbs(5000, s, seed=42) for s in surveys}

    # Plot redshift distributions
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))

    for ax, survey in zip(axes, surveys):
        df = dfs[survey]
        ax.hist(df['z'], bins=30, alpha=0.7, density=True)
        ax.set_xlabel('Redshift z')
        ax.set_ylabel('Density')
        ax.set_title(f'{survey}: median z = {df["z"].median():.2f}')

    plt.tight_layout()
    plt.savefig('frb_redshift_comparison.png')

    # Print summary statistics
    for survey in surveys:
        df = dfs[survey]
        print(f"\n{survey}:")
        print(f"  DM range: {df['DM'].min():.0f} - {df['DM'].max():.0f} pc/cm^3")
        print(f"  z range: {df['z'].min():.3f} - {df['z'].max():.3f}")
        print(f"  m_r range: {df['m_r'].min():.1f} - {df['m_r'].max():.1f}")
