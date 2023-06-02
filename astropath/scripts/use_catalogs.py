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

    import seaborn as sns
    from matplotlib import pyplot as plt

    from astropy import units

    from frb.surveys import survey_utils

    from astropath import path
    from astropath.scripts.utils import coord_arg_to_coord
    from astropath.utils import radec_to_coord

    scale = 0.5

    if pargs.ltype == 'ellipse':
        a, b, pa = [float(ip) for ip in pargs.lparam.split(',')]
        eellipse = {'a': a, 'b': b, 'theta': pa}

    # Load up the survey
    coord = radec_to_coord(coord_arg_to_coord(pargs.coord))
    survey = survey_utils.load_survey_by_name(
        pargs.survey, coord, pargs.ssize*units.arcmin)

    # Survey specific queries
    if pargs.survey == 'Pan-STARRS':
        query_fields = ['rPSFLikelihood']
    elif pargs.survey == 'DECaL':
        query_fields = ['shapedev_r', 'shapeexp_r']
    else:
        query_fields = None

    # Grab the catalo
    catalog = survey.get_catalog(query_fields=query_fields)

    if len(catalog) == 0:
        print(f"No objects in the catalog of your survey within {pargs.ssize} arcmin")
        return

    # Clean up the catalog
    if pargs.survey == 'Pan-STARRS':
        cut_size = catalog['rKronRad'] > 0.
        cut_mag = catalog['Pan-STARRS_r'] > 14. # Reconsider this
        cut_point = np.log10(np.abs(catalog['rPSFLikelihood'])) < (-2)
        keep = cut_size & cut_mag & cut_point
        # Half-light radius
        mag_key = 'Pan-STARRS_r'
        catalog['ang_size'] = catalog['rKronRad'].copy() 
    elif pargs.survey == 'DECaL':
        mag_key = 'DECaL_r'
        # Cuts
        cut_mag = (catalog[mag_key] > 14.) & np.isfinite(catalog[mag_key]) 
        cut_star = catalog['gaia_pointsource'] == 0
        keep = cut_mag & cut_star
        # Half-light radius
        catalog['ang_size'] = np.maximum(catalog['shapedev_r'], catalog['shapeexp_r'])
        zero = catalog['ang_size'] == 0.
        catalog['ang_size'][zero] = 1. # KLUDGE!!
    else:
        raise IOError(f"Not ready for this survey: {pargs.survey}")



    catalog = catalog[keep]

    if pargs.debug:
        if pargs.survey == 'Pan-STARRS':
            sns.histplot(x=catalog['Pan-STARRS_r'])#, bins=20)
            plt.show()
            sns.histplot(x=catalog['rKronRad'])#, bins=20)
            plt.show()
            sns.histplot(x=catalog['rPSFLikelihood'])#, bins=20)
            plt.show()
            embed(header='lowdm_bb: Need to set boxsize')

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
    Path.init_theta_prior('exp', 6., scale)

    # Priors
    p_O = Path.calc_priors()

    # Posterior
    P_Ox, P_Ux = Path.calc_posteriors('local', box_hwidth=box_hwidth, 
        max_radius=box_hwidth)

    # Print
    Path.candidates.sort_values(by='P_Ox', ascending=False, inplace=True)
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
