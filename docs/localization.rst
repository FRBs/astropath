**********************
Transient Localization
**********************

This doc describes the various ways one can define and
input the transient localization.

Data Model
==========

The localization is held internally in a *dict* 
with keys set by the data model defined in *localization.py*.

Error Ellipse
=============

The simplest approach is to define the localization as an
ellipse on the sky.  In this case, one inputs the center of
the ellipse, the semi-major axis ["a", in arcsec], 
the semi-minor axis ["b", in arcsec], and the PA on the sky (
"theta", defined in wacky astronomer fashion, i.e. deg E of N).

Here is an example::

    frb_coord = SkyCoord('21h44m25.255s -40d54m00.10s', frame='icrs')
    eellipse = dict(a=0.1, b=0.1, theta=0.)
    localiz = dict(type='eellipse', center_coord=frb_coord, eellipse=eellipse)
    assert localization.vette_localization(localiz)

The last line of code checks against the data model.

The code then defines the localization PDF as a 2D Gaussian
with sigma's "a" and "b" and orientation given by "theta".

Healpix
=======

For complex and/or large localizations, Healpix may offer the
best format.  Indeed, this is the preferred approach of the 
Graviational Wave community.

*astropath* accomodates two appraoches to defining the Healpix
localization.  We describe each in turn.

Nested
------

The first is termed "NESTED" and is a full Healpix
map of the sky with the PDF defined at every healpix pixel.

Here is an example using the healpix localization for
GW170817::

    lfile = os.path.join(resource_filename('astropath', 'data'), 'gw_examples',
                         'GW170817_skymap.fits.gz')
    gw170817 = hp.read_map(lfile)
    header = fits.open(lfile)[1].header
    #
    localiz = dict(type='healpix',
                   healpix_data=gw170817,
                   healpix_nside=header['NSIDE'],
                   healpix_ordering='NESTED',
                   healpix_coord='C')
    assert localization.vette_localization(localiz)


WCS
===