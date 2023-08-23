""" Methods related to public catalogs """

import numpy as np

from astropy import units

import seaborn as sns
from matplotlib import pyplot as plt

from frb.surveys import survey_utils

from IPython import embed

def query_catalog(psurvey:str, coord, ssize:float, debug:bool=False):
    """  Queary a public catalog for PATH analysis

    The catalog is cleaned up to remove objects that 
    are too bright or likely to be stars.

    Args:
        psurvey (str): Name of the survey ['Pan-STARRS', 'DECaL']
        coord (astropy.coordinates.SkyCoord): Coordinates of the transient
        ssize (float): Radius of the query in arcmin 
        debug (bool, optional): Show some plots. Defaults to False.

    Raises:
        IOError: _description_

    Returns:
        tuple: (astropy.table.Table) table of sources
            (str): Magnitude key (e.g. 'Pan-STARRS_r')
    """

    # Survey specific queries
    if psurvey == 'Pan-STARRS':
        query_fields = ['rPSFLikelihood']
    elif psurvey == 'DECaL':
        print("Using Pan-STARRS")
        query_fields = ['shapedev_r', 'shapeexp_r']
    else:
        query_fields = None
    print(f"Using: {psurvey} survey")

    survey = survey_utils.load_survey_by_name(
        psurvey, coord, ssize*units.arcmin)

    # Grab the catalo
    catalog = survey.get_catalog(query_fields=query_fields)

    # Empty?
    if len(catalog) == 0:
        print(f"No objects in the catalog of your survey within {ssize} arcmin")
        return catalog

    # Clean up the catalog
    if psurvey == 'Pan-STARRS':
        cut_size = catalog['rKronRad'] > 0.
        cut_mag = catalog['Pan-STARRS_r'] > 14. # Reconsider this
        cut_point = np.log10(np.abs(catalog['rPSFLikelihood'])) < (-2)
        keep = cut_size & cut_mag & cut_point
        # Half-light radius
        mag_key = 'Pan-STARRS_r'
        catalog['ang_size'] = catalog['rKronRad'].copy() 
        bad_ang = None
    elif psurvey == 'DECaL':
        mag_key = 'DECaL_r'
        # Cuts
        cut_mag = (catalog[mag_key] > 14.) & np.isfinite(catalog[mag_key]) 
        cut_star = catalog['gaia_pointsource'] == 0
        keep = cut_mag & cut_star
        # Half-light radius
        catalog['ang_size'] = np.maximum(catalog['shapedev_r'], catalog['shapeexp_r'])
        bad_ang = catalog['ang_size'] == 0.
        if np.any(bad_ang) > 0:
            print(f"WARNING:  Found {np.sum(bad_ang)} objects with zero ang_size. Setting to 1 arcsec")
        catalog['ang_size'][bad_ang] = 1.   # KLUDGE!!
    else:
        raise IOError(f"Not ready for this survey: {psurvey}")

    catalog = catalog[keep]

    if debug:
        if psurvey == 'Pan-STARRS':
            sns.histplot(x=catalog['Pan-STARRS_r'])#, bins=20)
            plt.show()
            sns.histplot(x=catalog['rKronRad'])#, bins=20)
            plt.show()
            sns.histplot(x=catalog['rPSFLikelihood'])#, bins=20)
            plt.show()
            embed(header='lowdm_bb: Need to set boxsize')

    return catalog, mag_key