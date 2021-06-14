""" PATH class to guide analysis """

import pandas
import numpy as np
import logging

from astropy.coordinates import SkyCoord

from astropath import candidates
from astropath import bayesian
from astropath import localization
from astropath import priors


class PATH(object):
    """Convenience class to run PATH analysis

    Attributes:
        candidates (pandas.DataFrame): 
            Table of host galaxy candidates
            Primary keys are:  
                ra: RA in deg
                dec: Dec in deg
                ang_size: Angular size of the galaxy in arcsec
                mag (optional): Apparent r-band magnitude
        cand_prior (dict):
            Holds the prior assumptions for the candidates
            see priors.py for further details
        theta_prior (dict):
            Holds the prior assumptions for transient distribution
            around the galaxy
            see priors.py for further details
        localize (dict):
            Holds the transient localization
            see localization.py for further details

    """
    def __init__(self):
        self.candidates = None  # pandas.DataFrame
        self.cand_prior = None  # dict
        self.theta_prior = None  # dict
        self.localiz = None  # dict

        # Extras
        self.raw_prior_Oi = None
        self.prior_Oi = None
        self.cand_coords = None
    
    def init_candidates(self, ra:np.ndarray, dec:np.ndarray, ang_size:np.ndarray, mag=None):
        """Generate the candidate table
        Also generate the self.cand_coords 

        Args:
            ra (np.ndarray): RA [deg]
            dec (np.ndarray): Dec [deg]
            ang_size (np.ndarray): Angular size [arcsec]
            mag ([type], optional): apparent magnitude [mag]
        """
        # Ingest
        self.candidates = pandas.DataFrame(dict(ra=ra, dec=dec, ang_size=ang_size))
        if mag is not None:
            self.candidates['mag'] = mag
        # Vet
        assert candidates.vet_candidates(self.candidates), 'Bad candidate input'
        # 
        self.init_cand_coords()
        logging.info("Candidates are ready!")

    def init_cand_coords(self):
        """Simple method to generate the coords object
        from the Table.  
        
        In case the candidates are set in some other fashion.
        """
        # Add SkyCoord objects
        self.cand_coords = SkyCoord(ra=self.candidates.ra, 
                                             dec=self.candidates.dec, 
                                             unit='deg')

    def init_cand_prior(self, P_O_method:str, P_U=0.):
        """Ingest the candidate prior approach and P_U

        Args:
            P_O_method (str): Method for calculation candidate prior
            P_U (float, optional): Unseen prior. Defaults to 0..

        Raises:
            IOError: [description]
        """
        if self.candidates is None:
            raise IOError("You must init candidates before their prior")
        self.cand_prior = dict(P_O_method=P_O_method, 
                               P_U=P_U)
        # Vet
        assert priors.vet_cand_prior(self.cand_prior, self.candidates), 'Bad candidate prior input'

    def init_theta_prior(self, PDF:str, max:float):
        """Ingest the theta (offset) prior

        Args:
            PDF (str): Method
            max (float): Maximum offset allowed
        """
        self.theta_prior = dict(PDF=PDF, max=max)
        # Vet
        assert priors.vet_theta_prior(self.theta_prior)
        logging.info("Priors are ready!")

    def init_localization(self, ltype:str, **kwargs):
        """Ingets the localization information

        Args:
            ltype (str): [description]
            kwargs: Other parameters definig the localization
                Depends on the ltype.  See localization.py
        """
        self.localiz = dict(type=ltype)
        self.localiz.update(kwargs)
        # Vet
        assert localization.vet_localization(self.localiz), 'Bad candidate prior input'
        logging.info("Localization is ready!")

    def calc_priors(self):
        """Calculate and normalize the P(O) values for the candidates

        Raises:
            IOError: [description]

        Returns:
            np.ndarray: P(O) values.  Also ingested in self.candidates['P_O']
        """
        if self.candidates is None or self.cand_prior is None:
            raise IOError("Init candidates and cand_prior first!")

        # Raw priors
        logging.info("Calculating priors")
        self.raw_prior_Oi = priors.raw_prior_Oi(
            self.cand_prior['P_O_method'], self.candidates['ang_size'], 
            mag=self.candidates['mag'] if 'mag' in self.candidates.keys() else None)

        # Normalize
        logging.info("Normalizing priors")
        self.prior_Oi = priors.renorm_priors(self.raw_prior_Oi, 
                                             self.cand_prior['P_U'])

        # Add to candidate table
        logging.info("Adding prior values to candidates table")
        self.candidates['P_O'] = self.prior_Oi

        # Return them too
        return self.prior_Oi

    def calc_posteriors(self, method:str, step_size=0.1, box_hwidth=None,
                        max_radius=None):
        """Calculate the posteriors

        Args:
            method (str): Approach
            step_size (float, optional): [description]. Defaults to 0.1.
            box_hwidth ([type], optional): [description]. Defaults to None.
            max_radius ([type], optional): [description]. Defaults to None.

        Raises:
            IOError: [description]
            ValueError: [description]
            IOError: [description]
            IOError: [description]

        Returns:
            [type]: [description]
        """
        # Check for inputs
        if self.candidates is None or self.cand_prior is None or self.localiz is None \
            or self.theta_prior is None: 
                raise IOError("Init everything first!!")
        # Check for priors!
        if 'P_O' not in self.candidates.keys():
            raise ValueError("You need to calculate the candidate priors first!!")

        # P(x|O)
        logging.info("Calculating p(x|O)")
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
            self.p_xOi = bayesian.px_Oi_local(
                        self.localiz, 
                        self.cand_coords,
                        self.candidates['ang_size'].values, 
                        self.theta_prior, 
                        step_size=step_size)
        self.candidates['p_xO'] = self.p_xOi

        # P(U|x)
        logging.info("Calculating p(x|U)")
        if self.cand_prior['P_U'] > 0.:
            if max_radius is None:
                raise IOError("Set max_radius given that P(U) > 0!!")
            self.p_xU = bayesian.px_U(max_radius)
        else:
            self.p_xU = 0.
        
        # p(x)
        logging.info("Calculating p(x)")
        self.p_x = self.cand_prior['P_U'] * self.p_xU + \
            np.sum(self.prior_Oi * self.p_xOi)

        # P(O|x)
        logging.info("Calculating P(O|x)")
        self.P_Oix = self.prior_Oi * self.p_xOi / self.p_x
        self.candidates['P_Ox'] = self.P_Oix

        # P(U|x)
        logging.info("Calculating P(U|x)")
        self.P_Ux = self.cand_prior['P_U'] * self.p_xU / self.p_x

        # Return
        return self.P_Oix, self.P_Ux
