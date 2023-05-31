""" Utility functions for astropath """

import numpy as np
import pandas
from pandas.core.frame import DataFrame


def match_ids(IDs, match_IDs, require_in_match=True):
    """ Match input IDs to another array of IDs (usually in a table)
    Return the rows aligned with input IDs
    Parameters
    ----------
    IDs : ndarray
    match_IDs : ndarray
    require_in_match : bool, optional
      Require that each of the input IDs occurs within the match_IDs
    Returns
    -------
    rows : ndarray
      Rows in match_IDs that match to IDs, aligned
      -1 if there is no match
    """
    rows = -1 * np.ones_like(IDs).astype(int)
    # Find which IDs are in match_IDs
    in_match = np.in1d(IDs, match_IDs)
    if require_in_match:
        if np.sum(~in_match) > 0:
            raise IOError("qcat.match_ids: One or more input IDs not in match_IDs")
    rows[~in_match] = -1
    #
    IDs_inmatch = IDs[in_match]
    # Find indices of input IDs in meta table -- first instance in meta only!
    xsorted = np.argsort(match_IDs)
    ypos = np.searchsorted(match_IDs, IDs_inmatch, sorter=xsorted)
    indices = xsorted[ypos]
    rows[in_match] = indices
    return rows


def vet_data_model(obj, dmodel:dict, verbose=True):
    """ Vet the input object against its data model

    Args:
        obj (dict or pandas.DataFrame):  Instance of the data model
        dmodel (dict): Data model
        verbose (bool): Print when something doesn't check

    Returns:
        tuple: chk (bool), disallowed_keys (list), badtype_keys (list)
    """

    chk = True
    # Loop on the keys
    disallowed_keys = []
    badtype_keys = []
    for key in obj.keys():
        # In data model?
        if not key in dmodel.keys():
            disallowed_keys.append(key)
            chk = False
            if verbose:
                print("Disallowed key: {}".format(key))

        # Check data type
        iobj = obj[key].values if isinstance(obj, pandas.DataFrame) else obj[key]
        if not isinstance(iobj,
                          dmodel[key]['dtype']):
            badtype_keys.append(key)
            chk = False        
            if verbose:
                print("Bad key type: {}".format(key))
    # Return
    return chk, disallowed_keys, badtype_keys



def radec_to_coord(radec, gal=False):
    """ Converts one of many of Celestial Coordinates
    `radec` formats to an astropy SkyCoord object. Assumes
    J2000 equinox.

    Parameters
    ----------
    radec : str or tuple or SkyCoord or list
        Examples:
        'J124511+144523',
        '124511+144523',
        'J12:45:11+14:45:23',
        ('12:45:11','+14:45:23')
        ('12 45 11', +14 45 23)
        ('12:45:11','14:45:23')  -- Assumes positive DEC
        (123.123, 12.1224) -- Assumed deg
        [(123.123, 12.1224), (125.123, 32.1224)]
    gal : bool, optional
      Input pair of floats are (l,b) in deg

    Returns
    -------
    coord : SkyCoord
      Converts to astropy.coordinate.SkyCoord (as needed)
      Returns a SkyCoord array if input is a list
    """
    from astropy.coordinates import SkyCoord
    if gal:
        frame = 'galactic'
    else:
        frame = 'icrs'

    # RA/DEC
    if isinstance(radec, (tuple)):
        if isinstance(radec[0], str):
            if radec[1][0] not in ['+', '-']:  #
                DEC = '+'+radec[1]
                warnings.warn("Assuming your DEC is +")
            else:
                DEC = radec[1]
            #
            coord = SkyCoord(radec[0]+DEC, frame=frame,
                                  unit=(u.hourangle, u.deg))
        else:
            if frame == 'galactic':
                coord = SkyCoord(l=radec[0], b=radec[1], frame=frame, unit='deg')
            else:
                coord = SkyCoord(ra=radec[0], dec=radec[1], frame=frame, unit='deg')
    elif isinstance(radec,SkyCoord):
        coord = radec
    elif isinstance(radec,str):
        # Find first instance of a number (i.e. strip J, SDSS, etc.)
        for ii in range(len(radec)):
            if radec[ii].isdigit():
                break
        radec = radec[ii:]
        #
        if ':' in radec:
            coord = SkyCoord(radec, frame='icrs', unit=(u.hourangle, u.deg))
        else:  # Add in :
            if ('+' in radec) or ('-' in radec):
                sign = max(radec.find('+'), radec.find('-'))
            else:
                raise ValueError("radec must include + or - for DEC")
            newradec = (radec[0:2]+':'+radec[2:4]+':'+radec[4:sign+3] +':'+radec[sign+3:sign+5]+':'+radec[sign+5:])
            coord = SkyCoord(newradec, frame='icrs', unit=(u.hourangle, u.deg))
    elif isinstance(radec,list):
        clist = []
        for item in radec:
            clist.append(radec_to_coord(item,gal=gal))
        # Convert to SkyCoord array
        ras = [ii.icrs.ra.value for ii in clist]
        decs = [ii.icrs.dec.value for ii in clist]
        return SkyCoord(ra=ras, dec=decs, unit='deg')
    else:
        raise IOError("Bad input type for radec")
    # Return
    return coord
