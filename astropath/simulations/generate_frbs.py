"""
Module for generating simulated FRB populations with properties (DM, z, Mr, mr)
for different surveys (CHIME, DSA, ASKAP).

This module samples from survey-specific P(DM,z) grids and host galaxy
magnitude distributions to generate realistic FRB populations.
"""

import numpy as np
import pandas
from scipy.interpolate import interp1d
from scipy import stats
import random

from astropy import units as u
from astropy.cosmology.realizations import Planck18

from frb.dm import prob_dmz
from frb.galaxies import hosts as hosts_mod

# Survey-specific telescope grid mappings
# Maps survey name to the grid filename used in frb.dm.prob_dmz
SURVEY_GRIDS = {
    'CHIME': 'CHIME_pzdm.npz',
    'DSA': 'DSA_pzdm.npz',
    'ASKAP': 'CRAFT_class_I_and_II_pzdm.npz',
    'CRAFT': 'CRAFT_class_I_and_II_pzdm.npz',
    'CRAFT_ICS_1300': 'CRAFT_ICS_1300_pzdm.npz',
    'CRAFT_ICS_892': 'CRAFT_ICS_892_pzdm.npz',
    'CRAFT_ICS_1632': 'CRAFT_ICS_1632_pzdm.npz',
    'Parkes': 'parkes_mb_class_I_and_II_pzdm.npz',
    'FAST': 'FAST_pzdm.npz',
}

# Default cosmology for distance modulus calculations
DEFAULT_COSMO = Planck18


def _build_cumulative_interpolator(values, pdf):
    """
    Build a cumulative distribution interpolator for inverse transform sampling.

    Args:
        values (np.ndarray): Array of values (e.g., DM or Mr)
        pdf (np.ndarray): Probability density values at each point

    Returns:
        scipy.interpolate.interp1d: Interpolator mapping uniform [0,1] -> values
    """
    cum = np.cumsum(pdf)
    cum[0] = 0.
    cum /= cum[-1]  # Normalize
    return interp1d(cum, values, bounds_error=False, fill_value=(values[0], values[-1]))


def _build_z_interpolators(pzdm, zvals, dmvals):
    """
    Build interpolators for P(z|DM) at each DM value.

    For each DM bin, creates an interpolator that samples z from the
    cumulative distribution P(z|DM).

    Args:
        pzdm (np.ndarray): 2D array of P(z,DM), shape (n_z, n_DM)
        zvals (np.ndarray): Redshift values
        dmvals (np.ndarray): DM values

    Returns:
        list: List of interpolators, one per DM bin
    """
    # Cumulative sum along z axis for each DM
    cum_all = np.cumsum(pzdm, axis=0)
    # Normalize each column
    norm = np.outer(np.ones(zvals.size), cum_all[-1, :])
    # Avoid division by zero
    norm[norm == 0] = 1.
    cum_all /= norm
    cum_all[0, :] = 0.

    # Build interpolators for each DM bin
    interpolators = []
    for ii in range(dmvals.size):
        # Handle edge case where all probabilities are zero
        if cum_all[-1, ii] == 0:
            interpolators.append(lambda x, z=zvals[0]: z)
        else:
            interpolators.append(
                interp1d(cum_all[:, ii], zvals,
                        bounds_error=False,
                        fill_value=(zvals[0], zvals[-1]))
            )
    return interpolators


def sample_dm_from_catalog(dm_values, n_samples, dm_range:tuple=None,
                           n_kde_points=500, seed=None):
    """
    Sample DM values from a KDE fit to observed catalog DMs.

    Args:
        dm_values (np.ndarray): Observed DM values from a catalog
        n_samples (int): Number of samples to generate
        dm_range (tuple): (min, max) DM range for KDE evaluation
        n_kde_points (int): Number of points for KDE evaluation
        seed (int, optional): Random seed for reproducibility

    Returns:
        np.ndarray: Sampled DM values
    """
    if dm_range is None:
        dm_range=(0., dm_values.max())
    if seed is not None:
        np.random.seed(seed)

    # Build KDE from observed DMs
    kernel = stats.gaussian_kde(dm_values)
    dm_grid = np.linspace(dm_range[0], dm_range[1], n_kde_points)
    dm_pdf = kernel(dm_grid)

    # Build interpolator and sample
    f_dm = _build_cumulative_interpolator(dm_grid, dm_pdf)
    rand = np.random.uniform(size=n_samples)
    return f_dm(rand)


def sample_redshifts_from_grid(dm_samples, pzdm, zvals, dmvals, seed=None):
    """
    Sample redshifts for given DM values using P(z|DM) grid.

    Args:
        dm_samples (np.ndarray): DM values to get redshifts for
        pzdm (np.ndarray): 2D probability grid P(z,DM), shape (n_z, n_DM)
        zvals (np.ndarray): Redshift values for the grid
        dmvals (np.ndarray): DM values for the grid
        seed (int, optional): Random seed for reproducibility

    Returns:
        np.ndarray: Sampled redshift values
    """
    if seed is not None:
        np.random.seed(seed)

    # Build z interpolators for each DM bin
    z_interpolators = _build_z_interpolators(pzdm, zvals, dmvals)

    # Sample redshifts
    n_samples = len(dm_samples)
    rand = np.random.uniform(size=n_samples)
    zs = np.zeros(n_samples)

    for kk, dm in enumerate(dm_samples):
        # Find closest DM bin
        imin = np.argmin(np.abs(dmvals - dm))
        zs[kk] = float(z_interpolators[imin](rand[kk]))

    return zs


def sample_host_Mr(n_samples, Mr_values=None, Mr_pdf=None,
                   Mr_range=(-25., -15.), n_kde_points=500, seed=None):
    """
    Sample host galaxy absolute r-band magnitudes.

    Can use either a provided (Mr, PDF) distribution or fit a KDE
    to provided Mr_values.

    Args:
        n_samples (int): Number of samples to generate
        Mr_values (np.ndarray, optional): Observed Mr values for KDE fitting
        Mr_pdf (tuple, optional): (Mr_array, pdf_array) pre-computed distribution
        Mr_range (tuple): (min, max) Mr range for KDE evaluation
        n_kde_points (int): Number of points for KDE evaluation
        seed (int, optional): Random seed for reproducibility

    Returns:
        np.ndarray: Sampled absolute magnitude values
    """
    if seed is not None:
        np.random.seed(seed)

    if Mr_pdf is not None:
        # Use provided PDF
        Mr_grid, pdf = Mr_pdf
    elif Mr_values is not None:
        # Build KDE from observed values
        kernel = stats.gaussian_kde(Mr_values)
        Mr_grid = np.linspace(Mr_range[0], Mr_range[1], n_kde_points)
        pdf = kernel(Mr_grid)
    else:
        raise ValueError("Must provide either Mr_values or Mr_pdf")

    # Build interpolator and sample
    f_Mr = _build_cumulative_interpolator(Mr_grid, pdf)
    rand = np.random.uniform(size=n_samples)
    return f_Mr(rand)


def calculate_apparent_mag(Mr, z, cosmo=None):
    """
    Calculate apparent magnitude from absolute magnitude and redshift.

    Args:
        Mr (np.ndarray): Absolute r-band magnitudes
        z (np.ndarray): Redshifts
        cosmo (astropy.cosmology, optional): Cosmology for distance modulus

    Returns:
        np.ndarray: Apparent r-band magnitudes
    """
    if cosmo is None:
        cosmo = DEFAULT_COSMO

    dist_mod = cosmo.distmod(z).value
    return dist_mod + Mr


def generate_frbs(n_frbs, survey, dm_catalog=None, 
    cosmo=None, seed=None, dm_range=None):
    """
    Generate a population of simulated FRBs with DM, z, Mr, and mr.

    This function generates FRBs by:
    1. Sampling DM from a KDE fit to observed catalog DMs (if provided)
       or directly from the P(DM,z) grid
    2. Sampling redshifts from survey-specific P(z|DM) grids
    3. Sampling host galaxy absolute magnitudes from known FRB host distribution
    4. Computing apparent magnitudes from z and Mr

    Args:
        n_frbs (int): Number of FRBs to generate
        survey (str): Survey name ('CHIME', 'DSA', 'ASKAP', etc.)
        dm_catalog (np.ndarray, optional): Observed DM values from catalog.
            If None, samples directly from the P(DM,z) grid.
        cosmo (astropy.cosmology, optional): Cosmology for calculations.
            Defaults to Planck18.
        seed (int, optional): Random seed for reproducibility

    Returns:
        pandas.DataFrame: DataFrame with columns:
            - 'DM': Extragalactic dispersion measure (pc/cm^3)
            - 'z': Redshift
            - 'M_r': Host galaxy absolute r-band magnitude
            - 'm_r': Host galaxy apparent r-band magnitude

    Raises:
        ValueError: If survey is not recognized

    Example:
        >>> df = generate_frbs(1000, 'CHIME')
        >>> print(df.head())
    """

    if survey not in SURVEY_GRIDS:
        raise ValueError(f"Unknown survey: {survey}. "
                        f"Available surveys: {list(SURVEY_GRIDS.keys())}")

    if cosmo is None:
        cosmo = DEFAULT_COSMO

    # Set master seed if provided
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    # Load the survey-specific P(z,DM) grid
    # First load CHIME grid to get z and DM arrays (they're the same for all)
    grid_data = prob_dmz.grab_repo_grid(SURVEY_GRIDS[survey])
    zvals = grid_data['z']
    dmvals = grid_data['DM']
    pzdm = grid_data['pzdm']

    # Step 1: Sample DM values
    print("Sampling DM values")
    if dm_catalog is not None:
        # Sample from KDE of observed DMs
        dm_samples = sample_dm_from_catalog(
            dm_catalog, n_frbs,
            dm_range=dm_range,
            seed=seed
        )
    else:
        # Sample directly from the P(DM,z) grid 
        grid_dict = {'pzdm': pzdm, 'z': zvals, 'DM': dmvals}
        df_temp = gen_random_FRBs(grid_dict, n_frbs)#, seed=seed)
        dm_samples = df_temp['DM'].values

    # Step 2: Sample redshifts given DMs
    print("Sampling redshifts")
    if dm_catalog is not None:
        # Need to sample z from P(z|DM) for each DM
        zs = sample_redshifts_from_grid(dm_samples, pzdm, zvals, dmvals)#, seed=seed)
    else:
        # Already sampled z along with DM
        zs = df_temp['z'].values

    # Step 3: Sample host galaxy absolute magnitudes
    # Load the Mr PDF from frb package
    print("Sampling host galaxy absolute magnitudes")
    Mr_grid, Mr_pdf_vals = hosts_mod.load_Mr_pdf()
    Mr_samples = sample_host_Mr(
        n_frbs,
        Mr_pdf=(Mr_grid, Mr_pdf_vals),
        #seed=seed
    )

    # Step 4: Calculate apparent magnitudes
    mr_samples = calculate_apparent_mag(Mr_samples, zs, cosmo=cosmo)

    # Build output DataFrame
    df_frbs = pandas.DataFrame({
        'DMeg': dm_samples,
        'z': zs,
        'M_r': Mr_samples,
        'm_r': mr_samples
    })

    return df_frbs

def gen_random_FRBs(grid:dict, nFRBs:int, seed:int=None):
    """
    Generate random Fast Radio Bursts (FRBs) based on a given probability grid.

    Parameters:
    -----------
    grid : dict
        A dictionary containing the probability grid with keys 'pzdm', 'z', and 'DM'.
        - 'pzdm' : 2D array-like, probability distribution over redshift (z) and dispersion measure (DM).
        - 'z' : 1D array-like, redshift values.
        - 'DM' : 1D array-like, dispersion measure values.
    nFRBs : int
        The number of random FRBs to generate.
    seed : int, optional
        Seed for the random number generator to ensure reproducibility. Default is None.

    Returns:
    --------
    pandas.DataFrame
        A DataFrame containing the generated FRBs with columns 'z' (redshift) and 'DM' (dispersion measure).
    """

    # Seed?
    if seed is not None:
        np.random.seed(seed)

    # Flatten 
    pzDM = grid['pzdm'].flatten()

    # Cum sum
    cum_sum = np.cumsum(pzDM)
    cum_sum /= cum_sum[-1]  # Normalize

    # Random numbers
    randu = np.random.uniform(size=nFRBs)

    # Assign to pzDM
    uidx = []
    for irand in randu:
        uidx.append(np.argmin(np.abs(irand-cum_sum)))
    # Unravel
    idx = np.unravel_index(uidx, grid['pzdm'].shape)

    # Generate the arrays
    z = grid['z'][idx[0]]
    DM = grid['DM'][idx[1]]

    # Pandas table
    df = pandas.DataFrame({'z':z, 'DM':DM})

    # Return
    return df