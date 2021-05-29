""" PATH class to guide analysis """

import pandas
import numpy as np

from astropy.coordinates import SkyCoord

from astropath import candidates
from astropath import utils
from astropath import bayesian
from astropath import localization
from astropath import priors


class PATH(object):
    def __init__(self):
        self.candidates = None  # pandas.DataFrame
        self.cand_prior = None  # dict
        self.theta_prior = None  # dict
        self.localiz = None  # dict

        # Extras
        self.raw_prior_Oi = None
        self.prior_Oi = None
        self.cand_coords = None
    
    def init_candidates(self, ra, dec, ang_size, mag=None):
        # Ingest
        self.candidates = pandas.DataFrame(dict(ra=ra, dec=dec, ang_size=ang_size))
        if mag is not None:
            self.candidates['mag'] = mag
        # Vet
        assert candidates.vet_candidates(self.candidates), 'Bad candidate input'
        # 
        self.init_cand_coords()

    def init_cand_coords(self):
        # Add SkyCoord objects
        self.cand_coords = SkyCoord(ra=self.candidates.ra, 
                                             dec=self.candidates.dec, 
                                             unit='deg')

    def init_cand_prior(self, P_O_method, P_U=0.):
        if self.candidates is None:
            raise IOError("You must init candidates before their prior")
        self.cand_prior = dict(P_O_method=P_O_method, 
                               P_U=P_U)
        # Vet
        assert priors.vet_cand_prior(self.cand_prior, self.candidates), 'Bad candidate prior input'

    def init_theta_prior(self, PDF, max):
        self.theta_prior = dict(PDF=PDF, max=max)
        # Vet
        assert priors.vet_theta_prior(self.theta_prior)

    def init_localization(self, ltype, **kwargs):
        self.localiz = dict(type=ltype)
        self.localiz.update(kwargs)
        # Vet
        assert localization.vet_localization(self.localiz), 'Bad candidate prior input'

    def calc_priors(self):
        if self.candidates is None or self.cand_prior is None:
            raise IOError("Init candidates and cand_prior first!")

        # Raw priors
        self.raw_prior_Oi = priors.raw_prior_Oi(
            self.cand_prior['P_O_method'], self.candidates['ang_size'], 
            mag=self.candidates['mag'] if 'mag' in self.candidates.keys() else None)

        # Normalize
        self.prior_Oi = priors.renorm_priors(self.raw_prior_Oi, 
                                             self.cand_prior['P_U'])

        # Add to candidate table
        self.candidates['P_O'] = self.prior_Oi

        # Return them too
        return self.prior_Oi

    def calc_posteriors(self, method, step_size=0.1, box_hwidth=None,
                        max_radius=None):
        # Check for inputs
        if self.candidates is None or self.cand_prior is None or self.localiz is None \
            or self.theta_prior is None: 
                raise IOError("Init everything first!!")
        # Check for priors!
        if 'P_O' not in self.candidates.keys():
            raise ValueError("You need to calculate the candidate priors first!!")

        # P(x|O)
        if method == 'fixed':
            if box_hwidth is None:
                raise IOError("Set box_hwidth to use method=fixed!")
            self.p_xOi = bayesian.px_Oi_fixedgrid(box_hwidth,  # box radius for grid, in arcsec
                        self.localiz, 
                        self.cand_coords,
                        self.candidates['ang_size'].values, 
                        self.theta_prior, 
                        step_size=step_size)
        elif method == 'local':
            pass
        self.candidates['p_xO'] = self.p_xOi

        # P(U|x)
        if self.cand_prior['P_U'] > 0.:
            if max_radius is None:
                raise IOError("Set max_radius given that P(U) > 0!!")
            self.p_xU = bayesian.px_U(max_radius)
        else:
            self.p_xU = 0.
        
        # p(x)
        self.p_x = self.cand_prior['P_U'] * self.p_xU + \
            np.sum(self.prior_Oi * self.p_xOi)

        # P(O|x)
        self.P_Oix = self.prior_Oi * self.p_xOi / self.p_x
        self.candidates['P_Ox'] = self.P_Oix

        # P(U|x)
        self.P_Ux = self.cand_prior['P_U'] * self.p_xU / self.p_x

        # Return
        return self.P_Oix, self.P_Ux
