""" Definitions and methods related to COSMOS """

import tqdm
import numpy as np

import pandas

from astropy import units
from astropy.coordinates import SkyCoord 

from IPython import embed

def cosmos_defs():
    cdefs = dict(plate_scale=0.05,
                 filter='ACS_i',
                 )
    return cdefs

def load_galaxies(cosmos_file:str='cosmos_acs_iphot_200709.feather'):
    # Load
    df = pandas.read_feather(cosmos_file)

    # Cut
    df = df[df['mu_class'] == 1]
    df = df[~np.isnan(df['kron_radius'])]

    pix_to_arcsec_kron = 0.03*units.arcsec
    df['half_light'] = df['kron_radius']*pix_to_arcsec_kron #arcsec

    # Return
    return df


def gen_frb_galaxy_pairs(cosmos_df:pandas.DataFrame,
    a_sigma:units.Quantity=1*units.arcsec,
    b_sigma:units.Quantity=1*units.arcsec,
    pa_mode:str='fixed',
    pa_value:float=0.):
# pix_to_arcsec = 0.05*units.arcsec

    frbs = []
    frb_with_gal = []

    for i, row in tqdm.tqdm(cosmos_df.iterrows()):
        coord = SkyCoord(ra=row['ra']*units.degree, dec=row['dec']*u.degree)

        # 
        loc_uncertainity = np.random.normal(scale=loc_sigma.value)*units.arcsec
        loc_pa = float(np.random.uniform(size=1, low=0., high=360.))
        coord_with_loc_uncertainity = coord.directional_offset_by(loc_pa*units.deg, loc_uncertainity)

    #     theta_max = pix_to_arcsec*row['flux_radius']
        theta_max = 2*row['half_light']
        galaxy_offset = np.random.uniform(low=0, high=theta_max)*units.arcsec
        gal_pa = float(np.random.uniform(size=1, low=0., high=360.))
        frb_coord = coord_with_loc_uncertainity.directional_offset_by(gal_pa*units.deg, galaxy_offset)

        row['frb_ra'] = frb_coord.ra.value
        row['frb_dec'] = frb_coord.dec.value
        frbs.append([row['name'], row['frb_ra'], row['frb_dec'], loc_sigma.value])
        
        frb_with_gal.append(row)
