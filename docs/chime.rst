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

This method requires that you have downloaded the COSMOS galaxy
feather file and placed into the `data/COSMOS` directory.  
See the `README.md` file there for the link.

Run MC
------

Now run the montecarlo.  Use the `run_mc()` method in
`calculations/CHIME/GBO/py/monte_carlo.py`, e.g.::

    run_mc('frb_monte_carlo_3x15.csv', 'PATH_3x15.csv')

This will produce a file `PATH_3x15.csv` with the results.
The method `parse_PATH()` in `calculations/CHIME/GBO/py/analysis.py`
can parse that file and produce a summary of the results.