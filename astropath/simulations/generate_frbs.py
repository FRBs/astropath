"""
Module for generating simulated FRB populations with properties (DM, z, Mr, mr)
for different surveys (CHIME, DSA, ASKAP).

This module samples from survey-specific P(DM,z) grids and host galaxy
magnitude distributions to generate realistic FRB populations.
"""

import numpy as np
import pandas
from importlib.resources import files
from scipy.interpolate import interp1d
from scipy.integrate import cumulative_trapezoid
from scipy import stats
import random

from astropy import units
from astropy.cosmology.realizations import Planck18

from frb.dm import prob_dmz

import warnings

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
    cum = cumulative_trapezoid(pdf, values, initial=0.)   # uses actual dx
    cum /= cum[-1]
    return interp1d(cum, values, bounds_error=False, fill_value=(values[0], values[-1]))


def kde_dm(data, c=60., npts=1400, dm_max=3200.):
    """KDE for DM_EG via y = log10(DM + c). Returns (DM_axis, pdf_in_DM)."""
    y = np.log10(data + c)
    kernel = stats.gaussian_kde(y)
    x = np.concatenate([np.linspace(0., 50., 200, endpoint=False),
                        np.geomspace(50., dm_max, npts)])
    pdf = kernel(np.log10(x + c)) / ((x + c) * np.log(10))   # Jacobian
    return x, pdf


def sample_dm_from_catalog(dm_values, n_samples, dm_range:tuple=None,
                           n_kde_points=2000, seed=None):
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
    if rng is None:
        rng = np.random.default_rng(seed)
    if dm_range is None:
        dm_range=(0., dm_values.max())

    # Build KDE from observed DMs
    dm_grid, dm_pdf = kde_dm(dm_values, c=60., npts=n_kde_points)

    # Build interpolator and sample
    f_dm = _build_cumulative_interpolator(dm_grid, dm_pdf)
    return f_dm(rng.uniform(size=n_samples))


def _cell_edges(centres):
    """
    Convert cell-centre coordinates to the n+1 cell edges.

    Handles non-uniform grids by placing interior edges at the midpoints and
    reflecting the first and last half-widths.

    Args:
        centres (np.ndarray): Monotonically increasing cell centres, length n.

    Returns:
        np.ndarray: Edges, length n+1.
    """
    c = np.asarray(centres, dtype=float)
    if c.ndim != 1 or c.size < 2:
        raise ValueError(f"centres must be 1-D with >=2 entries, got shape {c.shape}")
    if np.any(np.diff(c) <= 0):
        raise ValueError("centres must be strictly increasing")
    mid = 0.5 * (c[1:] + c[:-1])
    return np.concatenate([[c[0] - (mid[0] - c[0])], mid, [c[-1] + (c[-1] - mid[-1])]])


def build_z_cdf(pzdm, zvals, z_floor=0.0):
    """
    Build the per-DM cumulative distribution of z, evaluated at CELL EDGES.

    Args:
        pzdm (np.ndarray): Probability mass grid, shape (n_z, n_DM). Does not
            need to be normalized; each column is normalized independently.
        zvals (np.ndarray): Redshift cell centres, length n_z.
        z_floor (float): Lower bound clamped onto the first edge. The CHIME grid
            starts at z = 0.01 with dz = 0.01, so the first edge is 0.005; set
            `z_floor=0.01` if you would rather not sample below the grid's
            stated minimum. Default 0.0 (no clamping beyond non-negativity).

    Returns:
        tuple:
            - z_edges (np.ndarray): length n_z + 1
            - cdf (np.ndarray): shape (n_z + 1, n_DM), cdf[0, :] = 0,
              cdf[-1, :] = 1 for usable columns
            - empty (np.ndarray): bool mask, True for DM columns carrying no
              probability at all
    """
    pzdm = np.asarray(pzdm, dtype=float)
    zvals = np.asarray(zvals, dtype=float)
    if pzdm.ndim != 2:
        raise ValueError(f"pzdm must be 2-D (n_z, n_DM), got shape {pzdm.shape}")
    if pzdm.shape[0] != zvals.size:
        raise ValueError(
            f"pzdm has {pzdm.shape[0]} z rows but zvals has {zvals.size} entries. "
            "The grid may be transposed."
        )
    if np.any(pzdm < 0):
        raise ValueError("pzdm contains negative values")

    z_edges = _cell_edges(zvals)
    z_edges[0] = max(z_edges[0], z_floor, 0.0)

    # pzdm is MASS per cell -> cumsum gives the CDF at the right edge of each
    # cell. Prepending a zero gives the CDF at the left edge of the first cell,
    # so cdf and z_edges line up index for index.
    cdf = np.vstack([np.zeros((1, pzdm.shape[1])), np.cumsum(pzdm, axis=0)])

    total = cdf[-1, :].copy()
    empty = total <= 0
    total[empty] = 1.0                      # avoid divide-by-zero; masked below
    cdf /= total[None, :]

    if empty.any():
        warnings.warn(
            f"{empty.sum()} of {pzdm.shape[1]} DM columns carry no probability; "
            "draws against them return NaN.",
            stacklevel=2,
        )
    return z_edges, cdf, empty


def sample_redshifts_from_grid(dm_samples, pzdm, zvals, dmvals,
                               z_floor=0.0, rng=None, seed=None):
    """
    Sample one redshift per DM by inverting P(z | DM).

    Each DM is assigned to the grid cell that CONTAINS it, then z is drawn by
    inverse-CDF with linear interpolation across the containing z cell (i.e. a
    uniform density within the cell), so the result is continuous rather than
    snapped to grid points.

    Args:
        dm_samples (array): Extragalactic DM per FRB.
        pzdm (np.ndarray): Probability mass grid, shape (n_z, n_DM).
        zvals (np.ndarray): Redshift cell centres.
        dmvals (np.ndarray): DM cell centres.
        z_floor (float): See `build_z_cdf`.
        rng (np.random.Generator, optional): Preferred over `seed`.
        seed (int, optional): Used only if `rng` is None. Note this seeds a
            local Generator and does NOT touch the global numpy state, unlike
            the original `np.random.seed`.

    Returns:
        np.ndarray: Sampled redshifts, same length as `dm_samples`. Entries are
        NaN where the corresponding DM column carries no probability.

    Example:
        >>> d = np.load('CHIME_pzdm.npz')
        >>> z = sample_redshifts_from_grid(dm_cat, d['pzdm'], d['z'], d['DM'],
        ...                                seed=42)
    """
    if rng is None:
        rng = np.random.default_rng(seed)

    dm_samples = np.asarray(dm_samples, dtype=float)
    dmvals = np.asarray(dmvals, dtype=float)
    if pzdm.shape[1] != dmvals.size:
        raise ValueError(
            f"pzdm has {pzdm.shape[1]} DM columns but dmvals has {dmvals.size} entries."
        )
    if not np.all(np.isfinite(dm_samples)):
        raise ValueError("dm_samples contains non-finite values")

    z_edges, cdf, empty = build_z_cdf(pzdm, zvals, z_floor=z_floor)

    # Containing DM cell (identical to nearest-centre on a uniform grid, but
    # correct on a non-uniform one). Values outside the grid clamp to the ends.
    dm_edges = _cell_edges(dmvals)
    col = np.clip(np.searchsorted(dm_edges, dm_samples, side='right') - 1,
                  0, dmvals.size - 1)

    n_out = int(np.sum((dm_samples < dm_edges[0]) | (dm_samples > dm_edges[-1])))
    if n_out:
        warnings.warn(
            f"{n_out} DM values fall outside the grid range "
            f"[{dm_edges[0]:.1f}, {dm_edges[-1]:.1f}]; clamped to the edge columns.",
            stacklevel=2,
        )

    u = rng.uniform(size=dm_samples.size)
    out = np.full(dm_samples.size, np.nan)

    # Vectorize within each distinct DM column: the number of distinct columns
    # is at most n_DM regardless of how many FRBs are drawn.
    for j in np.unique(col):
        if empty[j]:
            continue
        sel = np.flatnonzero(col == j)
        c = cdf[:, j]
        i = np.clip(np.searchsorted(c, u[sel], side='right') - 1, 0, c.size - 2)
        width = c[i + 1] - c[i]
        # width == 0 only inside a flat stretch, which searchsorted cannot land
        # on with side='right'; guard anyway.
        frac = np.where(width > 0, (u[sel] - c[i]) / np.where(width > 0, width, 1.0), 0.0)
        out[sel] = z_edges[i] + frac * (z_edges[i + 1] - z_edges[i])

    return out


def sample_host_Mr(n_samples, Mr_pdf=None,
                   Mr_range=(-25., -15.), n_kde_points=500, rng=None, seed=None):
    """
    Sample host galaxy absolute r-band magnitudes.

    Can use either a provided (Mr, PDF) distribution or fit a KDE
    to provided Mr_values.

    Args:
        n_samples (int): Number of samples to generate
        Mr_pdf (tuple, optional): (Mr_array, pdf_array) pre-computed distribution
        Mr_range (tuple): (min, max) Mr range for KDE evaluation
        n_kde_points (int): Number of points for KDE evaluation
        seed (int, optional): Random seed for reproducibility

    Returns:
        np.ndarray: Sampled absolute magnitude values
    """
    if rng is None:
        rng = np.random.default_rng(seed)

    if Mr_pdf is not None:
        # Use provided PDF
        Mr_grid, pdf = Mr_pdf
    else:
        print("Using Lz values to sample host galaxy absolute magnitudes")
        # Load up Lz values
        host_file = files('astropath.data') / 'frb_surveys' / 'Lz_host_data.csv'
        df= pandas.read_csv(host_file)
        # Scale mrs with z's to find distribution of Mrs
        mrs = np.array(df['r-band'])
        zs_mrs = np.array(df['redshift'])

        # Get luminosity distance
        ds = Planck18.luminosity_distance(zs_mrs).to(units.parsec).value

        # Calculate absolute magnitudes
        Mr_values = mrs - 5. * np.log10(ds) + 5

        # Build KDE from observed values
        kernel = stats.gaussian_kde(Mr_values)
        Mr_grid = np.linspace(Mr_range[0], Mr_range[1], n_kde_points)
        pdf = kernel(Mr_grid)

    # Build interpolator and sample
    f_Mr = _build_cumulative_interpolator(Mr_grid, pdf)
    return f_Mr(rng.uniform(size=n_samples))


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
    seed_seq = np.random.SeedSequence(seed)
    rng_dm, rng_z, rng_mr = (np.random.default_rng(s) for s in seed_seq.spawn(3))
    if seed is not None:
        random.seed(seed)

    # Load the survey-specific P(z,DM) grid
    # First load CHIME grid to get z and DM arrays (they're the same for all)
    # grid_data = prob_dmz.grab_repo_grid(SURVEY_GRIDS[survey])
    pzdm_james23_fn = files('astropath.data') / 'simulations' / 'CHIME_progenitors_james23_case_bestfit.npz'
    grid_data = np.load(pzdm_james23_fn)
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
            rng=rng_dm,
        )
    else:
        # Sample directly from the P(DM,z) grid 
        grid_dict = {'pzdm': pzdm, 'z': zvals, 'DM': dmvals}
        df_temp = gen_random_FRBs(grid_dict, n_frbs, rng=rng_dm)#, seed=seed)
        dm_samples = df_temp['DM'].values

    # Step 2: Sample redshifts given DMs
    print("Sampling redshifts")
    if dm_catalog is not None:
        # Need to sample z from P(z|DM) for each DM
        zs = sample_redshifts_from_grid(dm_samples, pzdm, zvals, dmvals, rng=rng_z)#, seed=seed)
    else:
        # Already sampled z along with DM
        zs = df_temp['z'].values

    # Step 3: Sample host galaxy absolute magnitudes
    # Load the Mr PDF from frb package
    print("Sampling host galaxy absolute magnitudes")
    #Mr_grid, Mr_pdf_vals = hosts_mod.load_Mr_pdf()
    Mr_samples = sample_host_Mr(
        n_frbs,
        rng=rng_mr,
        #Mr_pdf=(Mr_grid, Mr_pdf_vals),
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

def gen_random_FRBs(grid:dict, nFRBs:int, seed:int=None, rng=None):
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
    if rng is None:
        rng = np.random.default_rng(seed)

    # Flatten 
    pzDM = grid['pzdm'].flatten()

    # Cum sum
    cum_sum = np.cumsum(pzDM)
    cum_sum /= cum_sum[-1]  # Normalize

    # Random numbers
    randu = rng.uniform(size=nFRBs)

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

def load_chime_cat1_DMeg():
    """ 
    Load the DM_extra galactic component from CHIME Catalog 1

    Returns:
    --------
    np.ndarray: DM_extra galactic component
    """
    # Load CHIME Catalog 1
    chime_cat_file = files('astropath.data') / 'frb_surveys' / 'chimefrbcat1.csv'
    df_dr1 = pandas.read_csv(chime_cat_file)

    # Make cut based on bonsai S/N
    cut_snr = 12.
    snr_cut = df_dr1['bonsai_snr'] > cut_snr
    df_dr1 = df_dr1[snr_cut].copy()

    # Subtract MW ISM component
    DMeg = np.nanmean([df_dr1['dm_exc_ne2001'].values,
        df_dr1['dm_exc_ymw16'].values], axis=0) #- 100. # Minus MW halo component

    # Return
    return DMeg
