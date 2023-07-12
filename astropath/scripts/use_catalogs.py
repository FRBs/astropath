#!/usr/bin/env python
"""
This script runs PATH on a input localization using public catalog data
"""
from IPython import embed

def parser(options=None):
    import argparse
    # Parse
    parser = argparse.ArgumentParser(description='Script to run PATH on a localization')
    parser.add_argument("coord", type=str, help="Central coordinates of the localization, e.g. J081240.7+320809 or 122.223,-23.2322 or 07:45:00.47,34:17:31.1")
    parser.add_argument("lparam", type=str, help="Localization parameters, e.g. 0.5,0.3,45 for ellipse which give semi-major and semi-minor axes and PA (in deg; E from N)")
    parser.add_argument("--ltype", type=str, default='ellipse', help="Localization type [ellipse] FUTURE: wcs, healpix")
    parser.add_argument("-U", "--PU", type=float, default=0., help="Prior on unseen galaxies")
    parser.add_argument("-s", "--survey", type=str, default='Pan-STARRS',
                        help="Public survey to use for the analysis ['Pan-STARRS', 'DECaL']")
    parser.add_argument("--scale", type=float, default=0.5, help="Scale for length in exponential prior")
    parser.add_argument("--ssize", type=float, default=5., help='Size of the survey in arcmin')
    parser.add_argument("--debug", default=False, action="store_true", help="debug?")
    parser.add_argument("-o", "--outfile", type=str, help="Name of the output file.  Should end in .csv")

    if options is None:
        pargs = parser.parse_args()
    else:
        pargs = parser.parse_args(options)
    return pargs


def main(pargs):
    """ Run
    """
    import numpy as np

    from astropath import path
    from astropath.scripts.utils import coord_arg_to_coord
    from astropath.utils import radec_to_coord
    from astropath import catalogs

    if pargs.ltype == 'ellipse':
        a, b, pa = [float(ip) for ip in pargs.lparam.split(',')]
        eellipse = {'a': a, 'b': b, 'theta': pa}

    # Load up the survey
    coord = radec_to_coord(coord_arg_to_coord(pargs.coord))

    # Grab the catalog
    catalog, mag_key = catalogs.query_catalog(pargs.survey, coord, pargs.ssize)


    # Set boxsize accoring to the largest galaxy (arcsec)
    box_hwidth = max(30., 10.*np.max(catalog['ang_size']))

    # Turn into a cndidate table
    Path = path.PATH()

   # Set up localization
    Path.init_localization('eellipse', 
                           center_coord=coord, 
                           eellipse=eellipse)
    # Coords
    Path.init_candidates(catalog['ra'],
                         catalog['dec'],
                         catalog['ang_size'],
                         mag=catalog[mag_key])

    # Candidate prior
    Path.init_cand_prior('inverse', P_U=pargs.PU)

    # Offset prior
    Path.init_theta_prior('exp', 6., pargs.scale)

    # Priors
    p_O = Path.calc_priors()

    # Posterior
    P_Ox, P_Ux = Path.calc_posteriors('local', box_hwidth=box_hwidth, 
        max_radius=box_hwidth)

    # Finish
    Path.candidates.sort_values(by='P_Ox', ascending=False, inplace=True)

    # Print
    print(Path.candidates[['ra', 'dec', 'ang_size', 'mag', 'P_O', 'P_Ox']])
    print(f"P_Ux = {Path.candidates['P_Ux'][0]}")

    # Save?
    if pargs.outfile is not None:
        Path.candidates.to_csv(pargs.outfile, index=False)
        print(f"Wrote: {pargs.outfile}")

    # Return
    return Path

# Test
# astropath_catalog 128.6800541558221,66.01075020181487 11.,11.,0. -U 0.2 --survey Pan-STARRS -o tst.csv

# You should receive:

#             ra        dec  ang_size        mag       P_O           P_Ox
#1    128.685113  66.007416  11.64190  16.223600  0.088637   9.094031e-01
#2    128.689961  66.009989   3.45769  20.777100  0.000560   6.501750e-03
#0    128.670817  66.010080   2.02068  22.538300  0.000110   1.440940e-03
#7    128.693225  66.000910   5.14814  18.084600  0.009618   7.935219e-04
#4    128.664696  66.009369   2.73424  22.483200  0.000115   3.902263e-04
