*******
Scripts
*******

This doc describes the scripts that are included with the package.

astropath_catalog
=================

This script ingests a localization, pulls a set of candidates
from a specified public survey (e.g. Pan-STARRS), and then
calculates the priors and posteriors for each candidate.


WARNING:  This script requires that the FRB repository be installed
including all of the survey requirements as described here:
https://github.com/FRBs/FRB/blob/main/docs/installing.rst 

Here is the current usage::

    usage: astropath_catalog [-h] [--ltype LTYPE] [-U PU] [-s SURVEY]
                            [--ssize SSIZE] [--debug] [-o OUTFILE]
                            coord lparam

    Script to run PATH on a localization

    positional arguments:
    coord                 Central coordinates of the localization, e.g.
                            J081240.7+320809 or 122.223,-23.2322 or
                            07:45:00.47,34:17:31.1
    lparam                Localization parameters, e.g. 0.5,0.3,45 for ellipse
                            which give semi-major and semi-minor axes and PA (in
                            deg; E from N)

    options:
    -h, --help            show this help message and exit
    --ltype LTYPE         Localization type [ellipse] FUTURE: wcs, healpix
    -U PU, --PU PU        Prior on unseen galaxies
    -s SURVEY, --survey SURVEY
                            Public survey to use for the analysis ['Pan-STARRS',
                            'DECaL']
    --scale SCALE         Scale for length in exponential prior
    --ssize SSIZE         Size of the survey in arcmin
    --debug               debug?
    -o OUTFILE, --outfile OUTFILE
                            Name of the output file. Should end in .csv

And here is an example::

    astropath_catalog 128.6800,66.01075 11.,11.,0. -U 0.2 --survey Pan-STARRS -o tst.csv

In this example, we have specified the coordinates in decimal degrees (ra, dec)
and the localization as an ellipse with semi-major and semi-minor axes of 11 arcseconds.
We have chosen the Pan-STARRS survey and an unseen prior of P_U=0.2. 
The output for this example looks like::

                ra        dec  ang_size        mag       P_O           P_Ox
    1    128.685113  66.007416  11.64190  16.223600  0.089033   9.098311e-01
    2    128.689961  66.009989   3.45769  20.777100  0.000563   6.504810e-03
    0    128.670817  66.010080   2.02068  22.538300  0.000110   1.441618e-03
    7    128.693225  66.000910   5.14814  18.084600  0.009661   7.938954e-04
    6    128.661478  66.003011   4.27073  18.915100  0.003839   3.687519e-04


    P_Ux = 0.0829160047730884

The columns are as follows::

    (ra) Right ascension of the candidate in decimal degrees,
    (dec) Declination of the candidate in decimal degrees,
    (ang_size) Angular size of the candidate in arcseconds,
    (mag) Magnitude of the candidate,
    (P_O) Prior on the candidate,
    (P_Ox) Posterior on the candidate.
