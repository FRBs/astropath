""" Analysis code for CHIME/GBO calculations """
import numpy as np

import pandas

from frb.surveys import catalog_utils

def parse_PATH(path_file:str, frb_file:str):

    # Load
    ftry = pandas.read_csv(path_file, index_col=0)
    ftry.iFRB = ftry.iFRB.values.astype(int)
    
    # COSMOS
    cosmos = pandas.read_feather('cosmos_acs_iphot_200709.feather')

    # FRBs
    frbs = pandas.read_csv(frb_file)

    # Add
    mt = catalog_utils.match_ids(frbs.gal_ID.values, cosmos.index)
    frbs['mag'] = cosmos.iloc[mt]['mag_best'].values

    # PATH
    POxs = [np.nan]*len(frbs)
    IDs = [-1]*len(frbs)

    for iFRB in np.unique(ftry.iFRB.values):
        best_cand = ftry[ftry.iFRB == iFRB].iloc[0] 
        # Save
        POxs[iFRB] = best_cand.P_Ox
        IDs[iFRB] = int(best_cand.gal_ID) 

    frbs['P_Ox'] = POxs
    frbs['PATH_ID'] = IDs

    return frbs