import os, sys
import numpy as np

from scipy.interpolate import interp1d

import pandas

from astropy import units

from frb.associate import frbassociate
from frb.associate import frbs

# Local
sys.path.append(os.path.abspath("../Analysis/py"))
import associate_defs

from IPython import embed

# Additional (unpublished) configs)

gdb_path = os.getenv('FRB_GDB')
base_config = frbs.base_config.copy()
updates = dict(
    name='FRB180301',
    image_file=os.path.join(gdb_path, 'Realfast', 'Bhandari2021', 'FRB180301_NOTr_ac.fits'),
    cut_size = 34.,
    filter = 'NOTr',
    ZP = 33.05,
    deblend=True,
    cand_bright=18.5,
    cand_separation=9*units.arcsec,
    plate_scale=0.2138*units.arcsec,
)
frb180301 = {**base_config, **updates}  # Use | in 3.9

def get_candidates():
    # Analyze without priors
    all_frbAs = []
    for ss, frb_name in enumerate(associate_defs.frb_list):
        config = getattr(frbs, frb_name.lower())
        config['skip_bayesian'] = True
        frbA = frbassociate.run_individual(config)
        all_frbAs.append(frbA)
    # df
    df = pandas.DataFrame(dict(frb=associate_defs.frb_list, frbA=all_frbAs))
    # Return
    return df


def run(prior, scale=1000, mag_lim=25.5, frb_pre=None):

    # Grab the candidates
    if frb_pre is None:
        frb_pre = get_candidates()

    # Init
    model_mags = []
    model_theta = []
    max_PMix = []

    # Bayesian
    for frbA in frb_pre.frbA:
        # Set priors
        frbA.calc_priors(prior['U'], method=prior['O'])
        prior['theta']['ang_size'] = frbA.candidates['half_light'].values
        frbA.set_theta_prior(prior['theta'])

        # Calculate p(O_i|x)
        frbA.calc_POx()

        # Record
        for cc, cand in frbA.candidates.iterrows():
            N = int(np.round(cand.P_Ox * scale))
            if N == 0:
                continue
            model_mags += [cand[frbA.filter]]*N
            model_theta += [(cand.separation/cand.half_light)]*N

        # Unseen
        N = int(np.round(frbA.P_Ux * scale))
        model_mags += [mag_lim] * N
        # Misc
        max_PMix.append(np.max(frbA.candidates.P_Ox))

    # Arrays
    model_mags = np.array(model_mags)
    model_theta = np.array(model_theta)
    max_PMix = np.array(max_PMix)

    # Return
    return frb_pre, model_mags, model_theta, max_PMix


def load_f_mL():
    # Grab m(L) table
    df = pandas.read_table('../Data/galLF_vs_z.txt', index_col=False)

    # Interpolate
    f_mL = interp1d(df.z, df['m_r(L*)'])

    # Return
    return f_mL
