""" Analysis code for CHIME/GBO calculations """

import pandas

from frb.surveys import catalog_utils

def parse_PATH():

    # Load
    ftry = pandas.read_csv('first_try.csv', index_col=0)
    ftry.iFRB = ftry.iFRB.values.astype(int)
    
    # COSMOS
    cosmos = pandas.read_feather('cosmos_acs_iphot_200709.feather')

    # FRBs
    frbs = pandas.read_csv('frb_monte_carlo.csv')

    # Add
    mt = catalog_utils.match_ids(frbs.gal_ID.values, cosmos.index)
    frbs['mag'] = cosmos.iloc[mt]['mag_best'].values