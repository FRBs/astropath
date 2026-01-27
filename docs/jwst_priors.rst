************
JWST Priors
************

This document describes the support for JWST imaging data in the
astropath prior calculations, specifically for the F200W filter.

Overview
========

The astropath package now supports galaxy number counts from JWST
imaging observations in addition to traditional optical r-band
observations. This enables proper host galaxy association for
transients detected in JWST fields.

The implementation uses galaxy number counts from Windhorst et al.
2024, which provide surface densities as a function of F200W
magnitude.

F200W Filter Support
====================

The F200W filter (JWST/NIRCam) is now supported as an option when
calculating raw priors for candidate host galaxies. This is implemented
through the ``filter`` parameter in the ``raw_prior_Oi`` function.

Supported Filters
-----------------

Currently, astropath supports the following filters:

* ``'r'`` -- r-band (optical), using Driver et al. 2016 number counts
* ``'F200W'`` -- JWST/NIRCam F200W (2.0 μm), using Windhorst et al. 2024 number counts

Usage
=====

Basic Example
-------------

To calculate priors for JWST F200W observations, specify the filter
parameter when calling ``raw_prior_Oi``::

    from astropath import priors
    import numpy as np

    # Example galaxy candidates
    mag_F200W = np.array([22.5, 23.0, 24.5])  # F200W magnitudes
    ang_size = np.array([1.0, 0.8, 0.5])      # Angular sizes in arcsec

    # Calculate raw priors using F200W filter
    raw_priors = priors.raw_prior_Oi(
        method='inverse',
        ang_size=ang_size,
        mag=mag_F200W,
        filter='F200W'
    )

Prior Calculation Methods
--------------------------

The same prior calculation methods used for r-band work with F200W:

* ``'inverse'`` -- :math:`P(O_i) \propto 1/\Sigma_m`
* ``'inverse_ang'`` -- :math:`P(O_i) \propto 1/(\Sigma_m \cdot \phi)`
* ``'inverse_ang2'`` -- :math:`P(O_i) \propto 1/(\Sigma_m \cdot \phi^2)`
* ``'identical'`` -- All galaxies equally weighted

where :math:`\Sigma_m` is the surface density of galaxies brighter
than magnitude :math:`m`, and :math:`\phi` is the angular size.

Complete Workflow
-----------------

Here is a complete example including prior normalization::

    from astropath import priors
    import numpy as np

    # Define candidate galaxies with F200W magnitudes
    candidates = {
        'mag': np.array([22.0, 23.5, 24.0]),
        'ang_size': np.array([2.0, 1.5, 1.0])
    }

    # Calculate raw priors with F200W filter
    raw_priors = priors.raw_prior_Oi(
        method='inverse',
        ang_size=candidates['ang_size'],
        mag=candidates['mag'],
        filter='F200W'
    )

    # Normalize priors (P_U is prior for unseen host)
    P_U = 0.0
    normalized_priors = priors.renorm_priors(raw_priors, P_U)

    print("Normalized priors:", normalized_priors)

Implementation Details
======================

Galaxy Number Counts
--------------------

The F200W number counts are derived from Windhorst et al. 2024 and
are stored as a data file in the astropath package. The
``windhorst_sigma`` function in the ``chance`` module reads these
number counts and interpolates them using a cubic spline to provide
smooth surface density estimates at any F200W magnitude.

The function signature is::

    def windhorst_sigma(mag):
        """
        Estimated incidence of galaxies per sq arcsec with F200W > mag
        using Windhorst et al. 2024 number counts.

        Args:
            mag: F200W band magnitude of galaxy

        Returns:
            Galaxy number density (per sq arcsec)
        """

Filter Selection
----------------

The ``raw_prior_Oi`` function automatically selects the appropriate
number count model based on the ``filter`` parameter:

* ``filter='r'`` uses ``chance.driver_sigma()``
* ``filter='F200W'`` uses ``chance.windhorst_sigma()``

If an unsupported filter is specified, a ``ValueError`` will be raised
with a list of supported filters.

Comparison with Optical Priors
===============================

When working with JWST data, it's important to note that:

1. **Magnitude System**: F200W is in the AB magnitude system, same as
   r-band, so magnitudes can be compared directly in terms of flux.

2. **Wavelength Regime**: F200W (2.0 μm) probes rest-frame optical/NIR
   at higher redshifts compared to optical r-band (0.6 μm).

3. **Depth**: JWST typically reaches fainter magnitudes than ground-based
   r-band imaging, allowing detection of fainter potential host galaxies.

4. **Number Counts**: The surface density :math:`\Sigma_m` differs between
   filters due to different wavelength selection effects and survey depths.

References
==========

* Windhorst et al. 2024 -- JWST F200W galaxy number counts
* Driver et al. 2016 -- Optical r-band galaxy number counts (original implementation)

See Also
========

* :doc:`offset_function` -- Offset function :math:`p(\omega|O)` descriptions
* :doc:`localization` -- Transient localization options
