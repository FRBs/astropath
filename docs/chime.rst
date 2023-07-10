*****
CHIME
*****

This doc describes aspects of *astropath* related
to CHIME

GBO Montecarlo
==============

We have added code to simulate PATH performance with a 
putative CHIME+GBO system.  This Montecarlo uses the
COSMOS field as the source of all host galaxies.

Generate FRBs
-------------


Use the method `generate_frbs()` in `calculations/CHIME/GBO/py/monte_carlo.py`
to generate a random set of FRBs, e.g.::

    generate_frbs('frb_monte_carlo_3x15.csv', 
                  radec_sigma=(3., 15.), debug=False, 
                  plots=False, nsample=10000)

Here, `radec_sigma` specifies the 1-sigma precision of the
loclization in RA and Dec, in arcsec.  
