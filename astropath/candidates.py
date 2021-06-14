""" Define data model for candidates """

import numpy as np
import pandas

from astropath import utils

candidates_dmodel = {
    'ra': dict(dtype=(np.ndarray),
                help='RA values in deg'),
    'dec': dict(dtype=(np.ndarray),
                help='Dec values in deg'),
    'ang_size': dict(dtype=(np.ndarray, list),
                help='Angular sizes for the candidates'),
    'mag': dict(dtype=(np.ndarray, list),
                help='Apparent magnitudes of the sources;  ' \
                    +'assumes r-band and corrected for Galactic extinction'),
}

def vet_candidates(candidates:pandas.DataFrame):
    """Vet the candidates Table
    i.e. confirm it corresponds to the DataModel

    Args:
        candidates (pandas.DataFrame): candidate table

    Raises:
        IOError: [description]

    Returns:
        bool: True if ok
    """
    chk, disallowed_keys, badtype_keys = utils.vet_data_model(candidates, candidates_dmodel)

    # Required keys
    required_keys = ['ra', 'dec', 'ang_size']
    for key in required_keys:
        if key not in candidates.keys():
            chk = False
    
    # Return
    return chk
    
