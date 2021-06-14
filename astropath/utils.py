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