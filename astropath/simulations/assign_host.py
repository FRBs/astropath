"""
Module for assigning simulated FRBs to host galaxies by apparent magnitude (m_r).

This module implements magnitude-based host assignment that matches FRBs to
galaxies from a catalog by their apparent r-band magnitudes, and then simulates
galactocentric offsets due to both the intrinsic FRB distribution and the
localization region.

Localization regions may be specified either as a single (a, b, PA) ellipse
applied to every FRB, or as a POPULATION of ellipses from which one is drawn at
random (with replacement) per FRB.  See `_resolve_localizations` and
`localizations_from_errors`.
"""

import os
import numpy as np
import random
import pandas as pd
import warnings
from pathlib import Path
from typing import Tuple, Optional, List, Union, Sequence

from scipy.spatial import cKDTree

from astropy import units
from astropy.coordinates import SkyCoord, match_coordinates_sky

from IPython import embed

# (a, b, PA) triple, an (N, 3) array, or a DataFrame with a/b/PA columns
LocalizationLike = Union[
    Tuple[float, float, float],
    Sequence[Tuple[float, float, float]],
    np.ndarray,
    pd.DataFrame,
]


def load_galaxy_catalog(catalog_fn: str = 'combined_HSC_DECaLs_HECATE_galaxies_hecatecut.parquet'):
    """
    Load possible host galaxy catalog with filename catalog_fn located
    in the directory indicated by the FRB_APATH environmental variable.

    Args:
        catalog_fn (str): filename of the catalog of possible hosts
    """
    frb_apath = os.environ.get('FRB_APATH')

    if frb_apath is not None:
        catalog_path = Path(frb_apath) / catalog_fn

        if catalog_path.exists():
            print(f"Loading real galaxy catalog from:")
            print(f"  {catalog_path}")
            galaxies = pd.read_parquet(catalog_path)

            return galaxies
        else:
            raise ValueError(f"Catalog not found at {catalog_path}")
    else:
        raise ValueError("FRB_APATH not set")


# ---------------------------------------------------------------------------
# Localization handling
# ---------------------------------------------------------------------------

def localizations_from_errors(
    ra_err: np.ndarray,
    dec_err: np.ndarray,
    pa: float = 90.0,
    dm: np.ndarray = None,
    min_err: float = 0.0,
    max_err: float = np.inf,
    drop_nonpositive: bool = True,
) -> np.ndarray:
    """
    Build a localization population array from per-FRB RA/Dec uncertainties.

    Intended for turning an observed catalog (e.g. BaseCat2) into the
    `localization` argument of `assign_frbs_to_hosts`.

    The convention matches `_apply_localization_error`: `a` is the offset scale
    along position angle `PA` (degrees East of North) and `b` is the scale along
    PA + 90.  With the default ``pa=90``, `a` runs East-West and `b` runs
    North-South, so `a` should be the RA uncertainty and `b` the Dec
    uncertainty.  Both are treated as 1-sigma Gaussian widths, truncated at
    3 sigma, exactly as for a scalar localization.

    Args:
        ra_err (array): RA uncertainty per FRB (arcsec), already sky-projected
            (i.e. including any cos(dec) factor).
        dec_err (array): Dec uncertainty per FRB (arcsec).
        pa (float): Position angle assigned to every ellipse (deg E of N).
            Default 90, pairing `a` with RA.
        dm (array, optional): Dispersion measure of each real FRB, on the SAME
            scale as the simulated `frb_df['DM']` -- i.e. EXTRAGALACTIC DM, with
            the Milky Way disk and halo contributions already subtracted. If
            given, the returned array has 4 columns (DM, a, b, PA) and enables
            joint DM-localization sampling in `assign_frbs_to_hosts`, which
            preserves the real correlation between burst DM and localization
            quality. If omitted, the returned array has 3 columns and
            localizations are drawn independently of DM.
        min_err (float): Discard rows with either axis below this (arcsec).
            Useful for stripping sentinels such as -99.99.
        max_err (float): Discard rows with either axis above this (arcsec).
        drop_nonpositive (bool): Discard non-finite or non-positive rows.

    Returns:
        np.ndarray: (N, 3) array of (a, b, PA), or (N, 4) array of
        (DM, a, b, PA) when `dm` is supplied. Ready to pass as `localization`.

    Example:
        >>> cat = pd.read_csv('basecat2_results.csv')
        >>> ra_e = pd.to_numeric(cat['RA Error (arcsec)'], errors='coerce')
        >>> de_e = pd.to_numeric(cat['Dec Error (arcsec)'], errors='coerce')
        >>> # independent sampling
        >>> locs = localizations_from_errors(ra_e, de_e)
        >>> # joint (DM, a, b, PA) sampling -- note DM must be EXTRAGALACTIC
        >>> dm_ex = pd.to_numeric(cat['Struct-max DM (pc cm-3)'], errors='coerce') - dm_mw
        >>> locs = localizations_from_errors(ra_e, de_e, dm=dm_ex)
        >>> df = assign_frbs_to_hosts(frbs, galaxies, localization=locs, seed=42)
    """
    a = np.asarray(ra_err, dtype=float)
    b = np.asarray(dec_err, dtype=float)

    if a.shape != b.shape:
        raise ValueError(f"ra_err and dec_err must have the same shape, got {a.shape} and {b.shape}")

    keep = np.ones(a.shape, dtype=bool)
    if drop_nonpositive:
        keep &= np.isfinite(a) & np.isfinite(b) & (a > 0) & (b > 0)
    keep &= (a >= min_err) & (b >= min_err)
    keep &= (a <= max_err) & (b <= max_err)

    dm_arr = None
    if dm is not None:
        dm_arr = np.asarray(dm, dtype=float)
        if dm_arr.shape != a.shape:
            raise ValueError(
                f"dm must have the same shape as ra_err/dec_err, got {dm_arr.shape} and {a.shape}"
            )
        keep &= np.isfinite(dm_arr) & (dm_arr > 0)

    n_drop = int((~keep).sum())
    if n_drop:
        print(f"localizations_from_errors: dropped {n_drop}/{len(a)} invalid rows")

    a, b = a[keep], b[keep]
    if len(a) == 0:
        raise ValueError("No valid localizations remain after filtering")

    pa_col = np.full(len(a), float(pa))
    if dm_arr is None:
        return np.column_stack([a, b, pa_col])

    dm_arr = dm_arr[keep]
    print(
        f"localizations_from_errors: DM column supplied "
        f"(median {np.median(dm_arr):.0f} pc/cm3) -> joint (DM, a, b, PA) sampling enabled"
    )
    return np.column_stack([dm_arr, a, b, pa_col])


def _check_localization_values(arr: np.ndarray):
    """Validate an (N, 3) (a, b, PA) array; raise on non-finite or non-positive axes."""
    if not np.all(np.isfinite(arr)):
        bad = int((~np.isfinite(arr)).any(axis=1).sum())
        raise ValueError(f"Localization contains non-finite values in {bad} row(s)")
    if np.any(arr[:, :2] <= 0):
        bad = int((arr[:, :2] <= 0).any(axis=1).sum())
        raise ValueError(
            f"Localization semi-axes must be positive; {bad} row(s) have a<=0 or b<=0. "
            "Check for sentinel values such as -99.99 (localizations_from_errors "
            "strips these for you)."
        )


def _resolve_localizations(
    localization: LocalizationLike,
    n_frbs: int,
    dm: np.ndarray = None,
    dm_neighbors: int = 25,
    dm_match_scale: str = 'log',
    debug: bool = False,
) -> np.ndarray:
    """
    Normalize the `localization` argument to a per-FRB (n_frbs, 3) array of
    (a, b, PA).

    Three input forms are accepted:

    1. A single (a, b, PA) triple -> broadcast unchanged to every FRB. This is
       the original behaviour and involves no randomness.

    2. A population of ellipses with shape (N, 3), or a DataFrame with columns
       'a', 'b', 'PA' -> one whole row drawn at random WITH REPLACEMENT per FRB,
       independently of any FRB property.

    3. A population with shape (N, 4) as (DM, a, b, PA), or a DataFrame that
       additionally has a 'DM' column -> JOINT (DM, a, b, PA) SAMPLING. For each
       simulated FRB, the `dm_neighbors` real FRBs closest in DM are found and
       one of them is drawn uniformly, so the simulated sample inherits the real
       correlation between burst DM and localization quality (brighter, better
       localized bursts sit preferentially at low DM). Requires `dm`.

    In all population cases whole rows are sampled, so correlations among a, b
    and PA are preserved. Draws use the global numpy RNG, so seeding via
    `np.random.seed` in the calling function makes the result reproducible.

    Args:
        localization: Ellipse or ellipse population (see above).
        n_frbs (int): Number of FRBs needing a localization.
        dm (array, optional): Extragalactic DM of each simulated FRB, length
            n_frbs. Required for form 3, ignored otherwise.
        dm_neighbors (int): Size of the DM neighbourhood sampled from in form 3.
            1 gives strict nearest-neighbour matching (tightest correlation, but
            heavy reuse of individual real FRBs); larger values loosen the
            correlation and increase diversity. Default 25.
        dm_match_scale (str): 'log' (default) or 'linear'. DM spans more than a
            decade, so matching in log DM gives neighbourhoods of roughly
            constant fractional width rather than ones that are far too wide at
            low DM and too narrow at high DM.
        debug (bool): Print extra diagnostics.

    Returns:
        np.ndarray: (n_frbs, 3) float array of (a, b, PA) in (arcsec, arcsec, deg).
        The DM column, if supplied, is used for matching and then dropped.

    Raises:
        ValueError: If the input cannot be interpreted, is empty, contains
            non-finite or non-positive axes, or a 4-column population is given
            without `dm`.
    """
    # --- DataFrame: pull DM/a/b/PA columns (case-insensitive) ---
    if isinstance(localization, pd.DataFrame):
        cols = {c.lower(): c for c in localization.columns}
        missing = [k for k in ('a', 'b', 'pa') if k not in cols]
        if missing:
            raise ValueError(
                f"localization DataFrame missing columns {missing}; "
                f"expected 'a', 'b', 'PA' (any case, plus optional 'DM'), "
                f"got {list(localization.columns)}"
            )
        want = ([cols['dm']] if 'dm' in cols else []) + [cols['a'], cols['b'], cols['pa']]
        arr = localization[want].to_numpy(dtype=float)
    else:
        arr = np.asarray(localization, dtype=float)

    # --- Single triple -> broadcast (original behaviour) ---
    if arr.ndim == 1:
        if arr.size != 3:
            raise ValueError(
                f"Scalar localization must be (a, b, PA) with 3 elements, got {arr.size}"
            )
        _check_localization_values(arr[None, :])
        print(
            f"Localization: fixed ellipse a={arr[0]:.2f}\" b={arr[1]:.2f}\" PA={arr[2]:.1f} deg "
            f"for all {n_frbs} FRBs"
        )
        return np.repeat(arr[None, :], n_frbs, axis=0)

    # --- Population ---
    if arr.ndim != 2 or arr.shape[1] not in (3, 4):
        raise ValueError(
            f"Localization population must have shape (N, 3) for (a, b, PA) or "
            f"(N, 4) for (DM, a, b, PA), got {arr.shape}"
        )
    if arr.shape[0] == 0:
        raise ValueError("Localization population is empty")

    joint = arr.shape[1] == 4
    pop_dm = arr[:, 0] if joint else None
    ellipses = arr[:, 1:] if joint else arr

    _check_localization_values(ellipses)

    if joint:
        if not np.all(np.isfinite(pop_dm)) or np.any(pop_dm <= 0):
            raise ValueError(
                "Localization DM column contains non-finite or non-positive values; "
                "strip sentinels before passing (localizations_from_errors does this)."
            )
        if dm is None:
            raise ValueError(
                "A 4-column (DM, a, b, PA) localization population was supplied, but no "
                "simulated DM array was passed. The FRB DataFrame needs a 'DM' column "
                "(extragalactic DM) for joint sampling. Pass a 3-column population "
                "instead to sample localizations independently of DM."
            )
        idx = _sample_by_dm(
            pop_dm, np.asarray(dm, dtype=float), dm_neighbors, dm_match_scale, debug=debug
        )
    else:
        if dm is not None and debug:
            print(
                "  (localization population has no DM column; sampling independently of DM)"
            )
        idx = np.random.randint(0, ellipses.shape[0], size=n_frbs)
        print(
            f"Localization: sampling {n_frbs} ellipses with replacement from a population "
            f"of {ellipses.shape[0]}, independently of FRB properties"
        )

    out = ellipses[idx]

    med_a, med_b = np.median(out[:, 0]), np.median(out[:, 1])
    area = np.pi * out[:, 0] * out[:, 1]
    print(
        f"  sampled median a={med_a:.1f}\"  b={med_b:.1f}\"  "
        f"area percentiles (sq arcsec) [10,50,90] = "
        f"{np.percentile(area, [10, 50, 90]).round(0)}"
    )
    if ellipses.shape[0] < n_frbs / 10:
        warnings.warn(
            f"Localization population ({ellipses.shape[0]}) is much smaller than the number "
            f"of FRBs ({n_frbs}); sampled ellipses will repeat heavily and the effective "
            f"variance of the localization distribution will be underestimated.",
            stacklevel=2,
        )

    return out


def _sample_by_dm(
    pop_dm: np.ndarray,
    frb_dm: np.ndarray,
    k: int,
    scale: str = 'log',
    debug: bool = False,
) -> np.ndarray:
    """
    Draw one real-FRB index per simulated FRB, matched in dispersion measure.

    For each simulated DM, the `k` nearest real FRBs in DM are located and one is
    chosen uniformly at random. This preserves the observed relationship between
    DM and localization quality while keeping the draw stochastic.

    Args:
        pop_dm (np.ndarray): DM of each real FRB in the localization population.
        frb_dm (np.ndarray): Extragalactic DM of each simulated FRB.
        k (int): Neighbourhood size; clipped to the population size.
        scale (str): 'log' or 'linear' matching space.
        debug (bool): Print the realized DM-matching residuals.

    Returns:
        np.ndarray: Integer indices into the localization population, one per
        simulated FRB.
    """
    if frb_dm.ndim != 1:
        raise ValueError(f"Simulated DM array must be 1-D, got shape {frb_dm.shape}")
    if not np.all(np.isfinite(frb_dm)) or np.any(frb_dm <= 0):
        n_bad = int((~np.isfinite(frb_dm)).sum() + (frb_dm <= 0).sum())
        raise ValueError(
            f"Simulated DM contains {n_bad} non-finite or non-positive values; "
            "joint DM-localization sampling needs a clean extragalactic DM for every FRB."
        )
    if scale not in ('log', 'linear'):
        raise ValueError(f"dm_match_scale must be 'log' or 'linear', got {scale!r}")

    k_eff = int(np.clip(k, 1, len(pop_dm)))
    if k_eff != k:
        warnings.warn(
            f"dm_neighbors={k} exceeds the localization population size "
            f"({len(pop_dm)}); using {k_eff}.",
            stacklevel=3,
        )

    xf = np.log10(frb_dm) if scale == 'log' else frb_dm
    xp = np.log10(pop_dm) if scale == 'log' else pop_dm

    tree = cKDTree(xp[:, None])
    _, nbr = tree.query(xf[:, None], k=k_eff)
    nbr = np.atleast_2d(nbr.T).T if k_eff > 1 else nbr.reshape(-1, 1)

    pick = np.random.randint(0, k_eff, size=len(xf))
    idx = nbr[np.arange(len(xf)), pick]

    # --- Diagnostics: how well does the real DM distribution cover the simulated one? ---
    lo, hi = pop_dm.min(), pop_dm.max()
    frac_out = float(np.mean((frb_dm < lo) | (frb_dm > hi)))
    print(
        f"Localization: joint (DM, a, b, PA) sampling -- {len(xf)} FRBs matched to a "
        f"population of {len(pop_dm)} real FRBs, {k_eff} nearest in {scale} DM"
    )
    print(
        f"  simulated DM median {np.median(frb_dm):.0f}, real DM median "
        f"{np.median(pop_dm):.0f} pc/cm3"
    )
    if frac_out > 0.01:
        print(
            f"  NOTE: {100 * frac_out:.1f}% of simulated FRBs fall outside the real DM range "
            f"[{lo:.0f}, {hi:.0f}]; these are matched to the nearest DM edge and their "
            f"localizations are extrapolated."
        )
    if debug:
        resid = np.abs(frb_dm - pop_dm[idx])
        print(
            f"  |DM_sim - DM_matched| percentiles [50,90,99] = "
            f"{np.percentile(resid, [50, 90, 99]).round(1)} pc/cm3"
        )
        n_used = len(np.unique(idx))
        print(f"  {n_used}/{len(pop_dm)} distinct real FRBs used")
        if np.std(pop_dm) > 0:
            r = np.corrcoef(frb_dm, np.pi * 1.0 * pop_dm[idx])[0, 1]
            print(f"  corr(DM_sim, DM_matched) = {r:.3f}")

    return idx


# ---------------------------------------------------------------------------
# Host assignment
# ---------------------------------------------------------------------------

def assign_frbs_to_hosts(
    frb_df: pd.DataFrame,
    galaxy_catalog: pd.DataFrame,
    localization: LocalizationLike,
    mag_range: Tuple[float, float] = None,  # (17., 28.),
    offset_function: str = 'exponential',
    scale: float = 0.5,
    trim_catalog: units.Quantity = 1 * units.arcmin,
    dm_neighbors: int = 25,
    dm_match_scale: str = 'log',
    seed: Optional[int] = None,
    debug: bool = False
) -> pd.DataFrame:
    """
    Assign FRBs to host galaxies based on apparent magnitude matching.

    This function assigns each FRB to a host galaxy by matching their apparent
    magnitudes (m_r). The matching algorithm uses a "fake coordinate" approach
    where magnitudes are encoded as declinations, allowing sky coordinate
    matching to effectively match by brightness.

    Each FRB is randomly placed within the host galaxy according to the indicated
    offset distribution (exponential by default), then the observed coordinates are offset
    according to the localization error ellipse.

    Args:
        frb_df (pd.DataFrame): FRB catalog with columns:
            - 'm_r': Apparent r-band magnitude of the host
            Additional columns are preserved in output
        galaxy_catalog (pd.DataFrame): Galaxy catalog with columns:
            - 'ra': Right ascension (degrees)
            - 'dec': Declination (degrees)
            - 'mag': Apparent r-band magnitude
            - 'half_light': Half-light radius (arcsec)
            - 'ID': Unique galaxy identifier
        localization: Error ellipse specification, in one of three forms:
            - a single (a, b, PA) triple applied to every FRB, where
              a = semi-major axis (arcsec), b = semi-minor axis (arcsec),
              PA = position angle (degrees East of North). Original behaviour,
              no randomness;
            - a POPULATION of ellipses -- an (N, 3) array, a sequence of
              triples, or a DataFrame with columns 'a', 'b', 'PA' -- from which
              one whole row is drawn at random WITH REPLACEMENT per FRB,
              independently of any FRB property;
            - a JOINT population -- an (N, 4) array of (DM, a, b, PA), or a
              DataFrame that also has a 'DM' column. Each simulated FRB is then
              matched to real FRBs of similar DM (see `dm_neighbors`) and given
              one of their localizations, so the simulated sample inherits the
              observed correlation between DM and localization quality. Requires
              a 'DM' column in `frb_df`, holding EXTRAGALACTIC DM on the same
              scale as the population's DM column.
            Use `localizations_from_errors` to build either population from an
            observed catalog's RA/Dec uncertainties.
            The realized per-FRB values are returned in the 'a', 'b' and 'PA'
            output columns.
        mag_range (tuple, optional): (min, max) magnitude range for FRB selection.
            FRBs outside this range are filtered out. Default: None (no cut)
        offset_function (str, optional): Function to use for generating galaxy positions (exponential, exponential_incorrect, uniform_1d, uniform_2d)
        scale (float, optional): Scale factor for exponential half-light radius or uniform distribution
            outer cutoff when offsetting FRBs due to intrinsic distribution.
            Smaller values concentrate FRBs closer to galaxy centers.
            Default: 0.5
        trim_catalog (units.Quantity, optional): Buffer to trim from catalog edges
            to ensure FRBs stay within analysis region. Default: 1 arcmin
        dm_neighbors (int, optional): Only used for joint (DM, a, b, PA)
            sampling. Number of real FRBs, nearest in DM, from which the
            localization is drawn uniformly. 1 is strict nearest-neighbour
            matching; larger values trade correlation strength for diversity.
            Default: 25
        dm_match_scale (str, optional): 'log' (default) or 'linear'; the space in
            which DM neighbours are found. Only used for joint sampling.
        seed (int, optional): Random seed for reproducibility. Also seeds the
            localization draw.
        debug (bool, optional): Enable debug output. Default: False

    Returns:
        pd.DataFrame: Table of FRB/host associations with columns:
            - 'ra': Observed FRB RA (degrees) - includes localization error
            - 'dec': Observed FRB Dec (degrees) - includes localization error
            - 'true_ra': True FRB RA in the galaxy (degrees)
            - 'true_dec': True FRB Dec in the galaxy (degrees)
            - 'gal_ID': ID of assigned host galaxy
            - 'gal_off': Offset from galaxy center (arcsec)
            - 'mag': Galaxy magnitude (m_r)
            - 'half_light': Galaxy half-light radius (arcsec)
            - 'loc_off': Localization error offset (arcsec)
            - 'FRB_ID': Original FRB index in input DataFrame
            - 'a': Localization semi-major axis for THIS FRB (arcsec)
            - 'b': Localization semi-minor axis for THIS FRB (arcsec)
            - 'PA': Localization position angle for THIS FRB (degrees)

    Raises:
        ValueError: If required columns are missing from input DataFrames

    Example:
        >>> frbs = generate_frbs(1000, 'CHIME')
        >>> galaxies = pd.read_parquet('galaxy_catalog.parquet')
        >>> # Fixed 25" circular localization
        >>> assignments = assign_frbs_to_hosts(frbs, galaxies, localization=(25, 25, 0))
        >>> # Realistic BaseCat2 localization distribution
        >>> locs = localizations_from_errors(cat['RA Error (arcsec)'],
        ...                                  cat['Dec Error (arcsec)'])
        >>> assignments = assign_frbs_to_hosts(frbs, galaxies, localization=locs)
        >>> # Joint DM-localization sampling (frbs must have a 'DM' column)
        >>> locs = localizations_from_errors(cat['RA Error (arcsec)'],
        ...                                  cat['Dec Error (arcsec)'], dm=dm_ex)
        >>> assignments = assign_frbs_to_hosts(frbs, galaxies, localization=locs)
    """
    # Set random seed if provided
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    # Validate input columns
    _validate_frb_columns(frb_df)
    _validate_galaxy_columns(galaxy_catalog)

    # Filter FRBs to reasonable magnitude range
    if mag_range is not None:
        mag_cut = (frb_df['m_r'] >= mag_range[0]) & (frb_df['m_r'] <= mag_range[1])
    else:
        mag_cut = np.ones(len(frb_df), dtype=bool)
    cut_frbs = frb_df[mag_cut].copy()

    if len(cut_frbs) == 0:
        raise ValueError(
            f"No FRBs remain after magnitude cut [{mag_range[0]}, {mag_range[1]}]. "
            f"Input m_r range: [{frb_df['m_r'].min():.2f}, {frb_df['m_r'].max():.2f}]"
        )

    print(f"Assigning {len(cut_frbs)} FRBs to hosts (filtered from {len(frb_df)})")

    # Draw a localization ellipse per FRB (or broadcast the single one given)
    loc_arr = _resolve_localizations(
        localization, len(cut_frbs),
        dm=_extract_dm(cut_frbs, localization),
        dm_neighbors=dm_neighbors, dm_match_scale=dm_match_scale, debug=debug,
    )

    # Trim catalog edges to maintain analysis region
    galaxy_cut = _trim_catalog(galaxy_catalog, trim_catalog)

    if len(galaxy_cut) == 0:
        raise ValueError("No galaxies remain after trimming catalog edges")

    # Match FRBs to galaxies by magnitude
    galaxy_indices = _match_by_magnitude(cut_frbs, galaxy_cut, debug=debug)
    galaxy_sample = galaxy_cut.loc[galaxy_indices]

    # Generate FRB positions within galaxies
    true_coords = _generate_galaxy_positions(
        galaxy_sample, scale=scale, function=offset_function,  # seed=seed
    )

    # Apply localization error
    obs_coords, loc_offsets = _apply_localization_error(
        true_coords, loc_arr,  # seed=seed
    )

    # Build output DataFrame
    df_out = _build_output_dataframe(
        obs_coords, true_coords, galaxy_sample,
        loc_offsets, loc_arr, cut_frbs.index.values
    )

    return df_out


def assign_frbs_random(
    frb_df: pd.DataFrame,
    galaxy_catalog: pd.DataFrame,
    localization: LocalizationLike,
    mag_range: Tuple[float, float] = None,
    offset_function: str = 'exponential',   # unused; kept for signature parity
    scale: float = 0.5,                     # unused; kept for signature parity
    trim_catalog: units.Quantity = 1 * units.arcmin,
    dm_neighbors: int = 25,
    dm_match_scale: str = 'log',
    seed: Optional[int] = None,
    debug: bool = False,
    coverage_catalog: pd.DataFrame = None,
    coverage_radius: units.Quantity = 0.5 * units.deg,
    include_boxes: List[Tuple[Tuple[float, float], Tuple[float, float]]] = None,
    max_batches: int = 200,
) -> pd.DataFrame:
    """
    Place FRB localization regions at random positions on-sky WITHIN THE FULL
    PRE-QUERIED FOOTPRINT, with NO association to a host galaxy.

    Null-hypothesis twin of `assign_frbs_to_hosts`, for measuring the rate of spurious
    high-confidence PATH associations (P(O|x) > 0.9) at random fields.

    The valid footprint is the UNION of:
        (1) disks of radius `coverage_radius` around the pre-queried Legacy Surveys /
            Pan-STARRS centers (~1 deg around bright galaxies), AND
        (2) any contiguous regions in `include_boxes` (e.g. the HSC-SSP XMM-LSS field).
    A localization center is accepted if it lands in EITHER (1) OR (2) -- not either/or
    at the call level, but a true combined footprint. Centers are drawn uniformly on the
    sphere (RA uniform, Dec uniform in sin(Dec)) over an enclosing box and rejection-
    tested against this union, so PATH always has catalog coverage over the analysis box.

    Args (differences from `assign_frbs_to_hosts`):
        localization: Single (a, b, PA) triple OR a population of ellipses
            sampled with replacement per FRB, exactly as in
            `assign_frbs_to_hosts`. Pass the SAME population you use for the
            host-assignment run so the null test carries the same localization
            distribution.
        coverage_catalog (pd.DataFrame, optional): Centers of the pre-queried disk
            patches (needs 'ra', 'dec'). Defaults to `galaxy_catalog`; pass the true
            query-center list if it differs, else the disk footprint may be undersized.
        coverage_radius (Quantity): Disk radius. Default 0.5 deg.
        include_boxes (list of ((ra_min, ra_max), (dec_min, dec_max)), optional):
            Contiguous fully-covered regions to ADD to the footprint, e.g.
            XMM-LSS: [((33., 38.), (-7., -2.))]. Boxes are edge-shrunk by `trim_catalog`.
            Boxes must not wrap RA=0 (a warning is issued if one appears to).
        offset_function, scale: ignored (no host to place the FRB within).

    Returns:
        pd.DataFrame with the SAME columns as `assign_frbs_to_hosts`:
            'ra','dec','a','b','PA','FRB_ID' -- POPULATED ('a','b','PA' now vary
                per FRB when a population is supplied)
            'true_ra','true_dec','gal_off','mag','half_light','loc_off' -- NaN
            'gal_ID' -- -99 int sentinel (build_digest tests this to skip host lookup)
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    _validate_frb_columns(frb_df)
    _validate_galaxy_columns(galaxy_catalog)

    if coverage_catalog is None:
        coverage_catalog = galaxy_catalog

    # Mirror the FRB subset assign_frbs_to_hosts would use (default None -> all FRBs)
    if mag_range is not None:
        mag_cut = (frb_df['m_r'] >= mag_range[0]) & (frb_df['m_r'] <= mag_range[1])
    else:
        mag_cut = np.ones(len(frb_df), dtype=bool)
    cut_frbs = frb_df[mag_cut].copy()

    if len(cut_frbs) == 0:
        raise ValueError(
            f"No FRBs remain after magnitude cut [{mag_range[0]}, {mag_range[1]}]. "
            f"Input m_r range: [{frb_df['m_r'].min():.2f}, {frb_df['m_r'].max():.2f}]"
        )
    n_needed = len(cut_frbs)
    print(f"Placing {n_needed} random localizations within pre-queried footprint "
          f"(filtered from {len(frb_df)})")

    # Draw a localization ellipse per FRB (or broadcast the single one given)
    loc_arr = _resolve_localizations(
        localization, n_needed,
        dm=_extract_dm(cut_frbs, localization),
        dm_neighbors=dm_neighbors, dm_match_scale=dm_match_scale, debug=debug,
    )

    # --- Disk geometry -------------------------------------------------------
    R = coverage_radius.to(units.deg).value
    buf = trim_catalog.to(units.deg).value
    R_eff = R - buf
    if R_eff <= 0:
        raise ValueError(
            f"coverage_radius ({R} deg) must exceed trim_catalog ({buf} deg)."
        )
    chord_thresh = 2.0 * np.sin(np.radians(R_eff) / 2.0)  # chord length for ang sep R_eff

    cov_xyz = _radec_to_unitvec(coverage_catalog['ra'].values,
                                coverage_catalog['dec'].values)
    tree = cKDTree(cov_xyz)

    # --- Box geometry (edge-shrunk by trim buffer) ---------------------------
    shrunk_boxes = []
    if include_boxes:
        for (ra_lo, ra_hi), (dec_lo, dec_hi) in include_boxes:
            if ra_hi - ra_lo > 180.:
                warnings.warn(
                    f"include_box RA span ({ra_lo},{ra_hi}) may wrap RA=0; "
                    "box membership test assumes no wrap."
                )
            # RA buffer scaled by cos(dec) so the shrink is a true angular buffer
            cosd = np.cos(np.radians(0.5 * (dec_lo + dec_hi)))
            ra_buf = buf / max(cosd, 1e-6)
            shrunk_boxes.append((
                (ra_lo + ra_buf, ra_hi - ra_buf),
                (dec_lo + buf, dec_hi - buf),
            ))

    # --- Enclosing sampling box (covers disks AND boxes) ---------------------
    # Full RA because the disk footprint wraps across RA=0; empty RA is rejected.
    ra_lo_s, ra_hi_s = 0.0, 360.0
    dec_lo_s = coverage_catalog['dec'].min() - R
    dec_hi_s = coverage_catalog['dec'].max() + R
    if shrunk_boxes:
        dec_lo_s = min(dec_lo_s, min(b[1][0] for b in shrunk_boxes))
        dec_hi_s = max(dec_hi_s, max(b[1][1] for b in shrunk_boxes))
    dec_lo_s = max(-90.0, dec_lo_s)
    dec_hi_s = min(90.0, dec_hi_s)
    sin_lo, sin_hi = np.sin(np.radians(dec_lo_s)), np.sin(np.radians(dec_hi_s))

    # --- Rejection sampling: uniform-on-sky within (disks UNION boxes) -------
    accepted_ra, accepted_dec = [], []
    n_acc, n_tried = 0, 0
    batch = max(4 * n_needed, 2000)
    for _ in range(max_batches):
        if n_acc >= n_needed:
            break
        cand_ra = np.random.uniform(ra_lo_s, ra_hi_s, size=batch)
        cand_dec = np.degrees(np.arcsin(np.random.uniform(sin_lo, sin_hi, size=batch)))
        dist, _ = tree.query(_radec_to_unitvec(cand_ra, cand_dec), k=1)
        keep = (dist <= chord_thresh) | _in_boxes(cand_ra, cand_dec, shrunk_boxes)
        accepted_ra.append(cand_ra[keep])
        accepted_dec.append(cand_dec[keep])
        n_acc += int(keep.sum())
        n_tried += batch
        p = max(n_acc / max(n_tried, 1), 1e-4)          # adapt to measured acceptance
        remaining = n_needed - n_acc
        batch = int(np.clip(1.5 * remaining / p, 2000, 5_000_000)) if remaining > 0 else batch
    else:
        raise RuntimeError(
            f"Only generated {n_acc}/{n_needed} localizations after {max_batches} batches "
            f"(acceptance ~{n_acc/max(n_tried,1):.4f}). Footprint may be tiny relative to "
            f"the sampling box, or coverage_radius/include_boxes too small."
        )

    rand_ra = np.concatenate(accepted_ra)[:n_needed]
    rand_dec = np.concatenate(accepted_dec)[:n_needed]

    if debug:
        p = n_acc / max(n_tried, 1)
        box_area = (ra_hi_s - ra_lo_s) * (sin_hi - sin_lo) * (180.0 / np.pi)  # sq deg
        # Fraction of accepted points that came via a box (overlap counts as box)
        if shrunk_boxes:
            in_box_final = _in_boxes(rand_ra, rand_dec, shrunk_boxes)
            print(f"{in_box_final.sum()}/{n_needed} localizations fell in include_boxes")
        print(f"acceptance ~{p:.4f}  ->  footprint area ~{p * box_area:.1f} sq deg")
        print(f"RA sampled [{ra_lo_s:.2f}, {ra_hi_s:.2f}], "
              f"Dec [{dec_lo_s:.3f}, {dec_hi_s:.3f}], R_eff={R_eff:.4f} deg")

    # --- Output (same columns / dummy scheme as assign_frbs_to_hosts) --------
    nan_col = np.full(n_needed, np.nan)
    return pd.DataFrame({
        'ra':         rand_ra,
        'dec':        rand_dec,
        'true_ra':    nan_col,
        'true_dec':   nan_col,
        'gal_ID':     np.full(n_needed, -99, dtype=int),   # sentinel: no host
        'gal_off':    nan_col,
        'mag':        nan_col,
        'half_light': nan_col,
        'loc_off':    nan_col,
        'FRB_ID':     cut_frbs.index.values,
        'a':          loc_arr[:, 0],
        'b':          loc_arr[:, 1],
        'PA':         loc_arr[:, 2],
    })


def assign_frbs_random_box(
    frb_df: pd.DataFrame,
    galaxy_catalog: pd.DataFrame,
    localization: LocalizationLike,
    mag_range: Tuple[float, float] = None,
    offset_function: str = 'exponential',   # unused; signature parity
    scale: float = 0.5,                     # unused; signature parity
    trim_catalog: units.Quantity = 1 * units.arcmin,
    dm_neighbors: int = 25,
    dm_match_scale: str = 'log',
    seed: Optional[int] = None,
    debug: bool = False,
    coverage_catalog: pd.DataFrame = None,
    coverage_box_width: units.Quantity = 0.7 * units.deg,
    include_boxes: List[Tuple[Tuple[float, float], Tuple[float, float]]] = None,
    max_batches: int = 200,
) -> pd.DataFrame:
    """
    Place FRB localization regions at random positions on-sky within the PATH
    footprint, with NO association to a host galaxy. Null-hypothesis twin of
    `assign_frbs_to_hosts` for measuring the spurious high-confidence rate
    (P(O|x) > 0.9).

    Footprint = (union of `coverage_box_width` x `coverage_box_width` boxes in RAW
    RA/Dec degrees around each `coverage_catalog` center, i.e. the HECATE hosts)
    UNION (`include_boxes`, e.g. the XMM-LSS field). Boxes are NOT cos(dec)-scaled,
    matching how the PATH catalog was queried. Centers are drawn uniformly on the
    sphere (RA uniform, Dec uniform in sin Dec) and accepted only inside the
    footprint, shrunk by `trim_catalog` so PATH always has coverage over the box.

    Args (differences from assign_frbs_to_hosts):
        localization: Single (a, b, PA) triple OR a population of ellipses
            sampled with replacement per FRB, exactly as in
            `assign_frbs_to_hosts`.
        coverage_catalog: centers of the per-host query boxes ('ra','dec'); the
            HECATE hosts. Defaults to galaxy_catalog. Pass the exact HECATE subset.
        coverage_box_width: full box width per host (default 0.7 deg).
        include_boxes: fully-covered regions to ADD, e.g. XMM-LSS
            [((33.575, 37.795), (-6.103, -3.191))].
        offset_function, scale: ignored (no host to place the FRB within).

    Returns:
        Same columns as assign_frbs_to_hosts:
            'ra','dec','a','b','PA','FRB_ID' -- POPULATED
            'true_ra','true_dec','gal_off','mag','half_light','loc_off' -- NaN
            'gal_ID' -- -99 int sentinel (build_digest tests this)
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    _validate_frb_columns(frb_df)
    _validate_galaxy_columns(galaxy_catalog)
    if coverage_catalog is None:
        coverage_catalog = galaxy_catalog

    if mag_range is not None:
        mag_cut = (frb_df['m_r'] >= mag_range[0]) & (frb_df['m_r'] <= mag_range[1])
    else:
        mag_cut = np.ones(len(frb_df), dtype=bool)
    cut_frbs = frb_df[mag_cut].copy()
    if len(cut_frbs) == 0:
        raise ValueError(
            f"No FRBs remain after magnitude cut [{mag_range[0]}, {mag_range[1]}]."
        )
    n_needed = len(cut_frbs)
    print(f"Placing {n_needed} random localizations (filtered from {len(frb_df)})")

    # Draw a localization ellipse per FRB (or broadcast the single one given)
    loc_arr = _resolve_localizations(
        localization, n_needed,
        dm=_extract_dm(cut_frbs, localization),
        dm_neighbors=dm_neighbors, dm_match_scale=dm_match_scale, debug=debug,
    )

    half = coverage_box_width.to(units.deg).value / 2.0
    buf = trim_catalog.to(units.deg).value
    half_eff = half - buf
    if half_eff <= 0:
        raise ValueError(
            f"trim_catalog ({buf} deg) must be < box half-width ({half} deg)."
        )

    # Chebyshev (box) membership: nearest center in (ra,dec) coord space, p=inf
    centers = np.column_stack([coverage_catalog['ra'].values,
                               coverage_catalog['dec'].values])
    tree = cKDTree(centers)

    shrunk_boxes = []
    if include_boxes:
        for (ra_lo, ra_hi), (dec_lo, dec_hi) in include_boxes:
            shrunk_boxes.append(((ra_lo + buf, ra_hi - buf),
                                 (dec_lo + buf, dec_hi - buf)))

    ra_lo_s, ra_hi_s = 0.0, 360.0
    dec_lo_s = coverage_catalog['dec'].min() - half
    dec_hi_s = coverage_catalog['dec'].max() + half
    if shrunk_boxes:
        dec_lo_s = min(dec_lo_s, min(b[1][0] for b in shrunk_boxes))
        dec_hi_s = max(dec_hi_s, max(b[1][1] for b in shrunk_boxes))
    dec_lo_s = max(-90.0, dec_lo_s)
    dec_hi_s = min(90.0, dec_hi_s)
    sin_lo, sin_hi = np.sin(np.radians(dec_lo_s)), np.sin(np.radians(dec_hi_s))

    accepted_ra, accepted_dec = [], []
    n_acc, n_tried = 0, 0
    batch = max(8 * n_needed, 5000)
    for _ in range(max_batches):
        if n_acc >= n_needed:
            break
        cand_ra = np.random.uniform(ra_lo_s, ra_hi_s, size=batch)
        cand_dec = np.degrees(np.arcsin(
            np.random.uniform(sin_lo, sin_hi, size=batch)))
        dist, _ = tree.query(np.column_stack([cand_ra, cand_dec]), k=1, p=np.inf)
        keep = (dist <= half_eff) | _in_boxes(cand_ra, cand_dec, shrunk_boxes)
        accepted_ra.append(cand_ra[keep])
        accepted_dec.append(cand_dec[keep])
        n_acc += int(keep.sum())
        n_tried += batch
        p = max(n_acc / max(n_tried, 1), 1e-5)
        remaining = n_needed - n_acc
        batch = int(np.clip(1.5 * remaining / p, 5000, 20_000_000)) if remaining > 0 else batch
    else:
        raise RuntimeError(
            f"Only generated {n_acc}/{n_needed} localizations after {max_batches} "
            f"batches (acceptance ~{n_acc/max(n_tried,1):.5f}). Footprint tiny vs box, "
            f"or coverage_box_width/include_boxes too small."
        )

    rand_ra = np.concatenate(accepted_ra)[:n_needed]
    rand_dec = np.concatenate(accepted_dec)[:n_needed]

    if debug:
        p = n_acc / max(n_tried, 1)
        box_area = (ra_hi_s - ra_lo_s) * (sin_hi - sin_lo) * (180.0 / np.pi)
        print(f"acceptance ~{p:.5f} -> footprint ~{p * box_area:.1f} sq deg")
        if shrunk_boxes:
            ib = _in_boxes(rand_ra, rand_dec, shrunk_boxes)
            print(f"{ib.sum()}/{n_needed} localizations in include_boxes")

    nan_col = np.full(n_needed, np.nan)
    return pd.DataFrame({
        'ra':         rand_ra,
        'dec':        rand_dec,
        'true_ra':    nan_col,
        'true_dec':   nan_col,
        'gal_ID':     np.full(n_needed, -99, dtype=int),
        'gal_off':    nan_col,
        'mag':        nan_col,
        'half_light': nan_col,
        'loc_off':    nan_col,
        'FRB_ID':     cut_frbs.index.values,
        'a':          loc_arr[:, 0],
        'b':          loc_arr[:, 1],
        'PA':         loc_arr[:, 2],
    })


def _radec_to_unitvec(ra_deg, dec_deg):
    """RA/Dec (deg) -> 3D unit vectors, shape (N, 3). Wrap-safe by construction."""
    ra = np.radians(np.asarray(ra_deg, dtype=float))
    dec = np.radians(np.asarray(dec_deg, dtype=float))
    cos_dec = np.cos(dec)
    return np.column_stack([cos_dec * np.cos(ra),
                            cos_dec * np.sin(ra),
                            np.sin(dec)])


def _in_boxes(ra_deg, dec_deg, boxes):
    """Boolean mask: True where (ra, dec) falls in ANY box. Boxes already edge-shrunk."""
    ra = np.asarray(ra_deg, dtype=float)
    dec = np.asarray(dec_deg, dtype=float)
    inside = np.zeros(ra.shape, dtype=bool)
    for (ra_lo, ra_hi), (dec_lo, dec_hi) in boxes:
        inside |= (ra >= ra_lo) & (ra <= ra_hi) & (dec >= dec_lo) & (dec <= dec_hi)
    return inside


def _needs_dm(localization: LocalizationLike) -> bool:
    """True if `localization` is a joint (DM, a, b, PA) population needing simulated DMs."""
    if isinstance(localization, pd.DataFrame):
        return 'dm' in {c.lower() for c in localization.columns}
    arr = np.asarray(localization, dtype=float)
    return arr.ndim == 2 and arr.shape[1] == 4


def _extract_dm(frb_df: pd.DataFrame, localization: LocalizationLike):
    """
    Pull the extragalactic DM column from the FRB table, but only when the
    supplied localization actually calls for joint sampling.

    Returns None for scalar or 3-column localizations, so the DM column is
    ignored unless it is needed.

    Raises:
        ValueError: If joint sampling is requested but `frb_df` has no 'DM'
            column.
    """
    if not _needs_dm(localization):
        return None
    if 'DMeg' not in frb_df.columns:
        raise ValueError(
            "Joint (DM, a, b, PA) localization sampling requires a 'DM' column in the FRB "
            f"DataFrame (extragalactic DM), but columns are {list(frb_df.columns)}. "
            "Pass a 3-column (a, b, PA) population to sample localizations independently."
        )
    return frb_df['DMeg'].values


def _validate_frb_columns(frb_df: pd.DataFrame):
    """Validate that FRB DataFrame has required columns."""
    required_cols = ['m_r']
    missing = [col for col in required_cols if col not in frb_df.columns]
    if missing:
        raise ValueError(f"FRB DataFrame missing required columns: {missing}")


def _validate_galaxy_columns(galaxy_df: pd.DataFrame):
    """Validate that galaxy DataFrame has required columns."""
    required_cols = ['ra', 'dec', 'mag', 'half_light', 'ID']
    missing = [col for col in required_cols if col not in galaxy_df.columns]
    if missing:
        raise ValueError(f"Galaxy catalog missing required columns: {missing}")


def _trim_catalog(galaxy_df: pd.DataFrame, trim: units.Quantity) -> pd.DataFrame:
    """
    Trim edges off catalog to maintain PATH analysis region.

    NOTE: this trims only the GLOBAL RA/Dec bounding box of the catalog, so for a
    footprint made of many separated islands it removes almost nothing and does
    not protect the per-island edges.

    Args:
        galaxy_df: Galaxy catalog
        trim: Buffer size to remove from edges

    Returns:
        Trimmed galaxy catalog
    """
    ra = galaxy_df.ra.values * units.deg
    dec = galaxy_df.dec.values * units.deg

    ra_min, ra_max = ra.min(), ra.max()
    dec_min, dec_max = dec.min(), dec.max()

    cut_ra = (ra > (ra_min + trim)) & (ra < (ra_max - trim))
    cut_dec = (dec > (dec_min + trim)) & (dec < (dec_max - trim))

    return galaxy_df[cut_ra & cut_dec]


def _match_by_magnitude(
    frb_df: pd.DataFrame,
    galaxy_df: pd.DataFrame,
    debug: bool = False
) -> np.ndarray:
    """
    Match FRBs to galaxies by apparent magnitude using fake coordinates.

    This implements a magnitude-matching algorithm from that creates
    "fake" sky coordinates where the declination encodes the magnitude,
    then uses astropy's match_coordinates_sky to match by brightness.

    The algorithm iteratively matches FRBs to galaxies, ensuring each galaxy
    is used only once.

    Args:
        frb_df: FRB catalog with 'm_r' column
        galaxy_df: Galaxy catalog with 'mag' column
        debug: Enable debug output

    Returns:
        Array of galaxy DataFrame indices matching each FRB
    """
    n_frbs = len(frb_df)

    # Create fake coordinates with magnitude encoded as declination
    # RA is set to 1 for all (doesn't matter, only dec is used for matching)
    fake_frb_coords = SkyCoord(
        ra=np.ones(n_frbs),
        dec=frb_df['m_r'].values,
        unit='deg'
    )

    fake_galaxy_coords = SkyCoord(
        ra=np.ones(len(galaxy_df)),
        dec=galaxy_df['mag'].values,
        unit='deg'
    )

    # Prepare for iterative matching
    galaxy_used = np.zeros(len(galaxy_df), dtype=bool)
    galaxy_indices = np.arange(len(galaxy_df))
    galaxy_df_indices = galaxy_df.index.values.copy()
    frb_assignments = -1 * np.ones(n_frbs, dtype=int)

    # Iteratively match FRBs to galaxies
    iteration = 0
    while np.any(frb_assignments < 0):
        iteration += 1
        n_remaining = np.sum(frb_assignments < 0)

        if debug or (iteration == 1) or (n_remaining < 100) or (iteration % 10 == 0):
            print(f"Iteration {iteration}: {n_remaining} FRBs remaining")
            print(f"  Brightest unassigned FRB: m_r = {np.min(fake_frb_coords[frb_assignments < 0].dec):.2f}")

        # Get unassigned FRBs and available galaxies
        unassigned_mask = frb_assignments < 0
        available_mask = ~galaxy_used

        sub_frb_coords = fake_frb_coords[unassigned_mask]
        sub_frb_indices = np.where(unassigned_mask)[0]

        sub_galaxy_coords = fake_galaxy_coords[available_mask]
        sub_galaxy_df_indices = galaxy_df_indices[available_mask]
        sub_galaxy_flag_indices = galaxy_indices[available_mask]

        # Check if we've run out of bright galaxies
        if np.max(sub_frb_coords.dec.deg) < np.min(sub_galaxy_coords.dec.deg):
            # Assign remaining FRBs to remaining galaxies by sorted magnitude
            print(f"Ran out of bright galaxies at iteration {iteration}")
            print(f"  Brightest remaining galaxy: m_r = {np.min(sub_galaxy_coords.dec.deg):.2f}")
            print(f"  Faintest remaining FRB: m_r = {np.max(sub_frb_coords.dec.deg):.2f}")

            srt_galaxies = np.argsort(sub_galaxy_coords.dec.deg)
            srt_frbs = np.argsort(sub_frb_coords.dec.deg)

            n_to_assign = min(len(srt_frbs), len(srt_galaxies))
            frb_assignments[sub_frb_indices[srt_frbs[:n_to_assign]]] = \
                sub_galaxy_df_indices[srt_galaxies[:n_to_assign]]

            if len(srt_frbs) > len(srt_galaxies):
                print(f"WARNING: {len(srt_frbs) - len(srt_galaxies)} FRBs could not be assigned")

            break

        # Match coordinates (effectively matching by magnitude)
        idx, d2d, _ = match_coordinates_sky(
            sub_frb_coords, sub_galaxy_coords, nthneighbor=1
        )

        if debug or iteration == 1:
            print(f"  Max magnitude separation: {d2d.max():.4f} deg")

        # Handle case where multiple FRBs match to same galaxy
        # Keep only first match for each unique galaxy
        unique_galaxies, unique_indices = np.unique(idx, return_index=True)

        # Assign these FRBs to their matched galaxies
        frb_assignments[sub_frb_indices[unique_indices]] = \
            sub_galaxy_df_indices[unique_galaxies]

        # Mark these galaxies as used
        galaxy_used[sub_galaxy_flag_indices[unique_galaxies]] = True

    print(f"Assignment complete after {iteration} iterations")

    # Verify all FRBs were assigned
    if np.any(frb_assignments < 0):
        n_unassigned = np.sum(frb_assignments < 0)
        raise RuntimeError(
            f"{n_unassigned} FRBs could not be assigned to galaxies. "
            "Galaxy catalog may be too small or magnitude distribution mismatch."
        )

    return frb_assignments


def _generate_galaxy_positions(
    galaxy_sample: pd.DataFrame,
    scale: float = 0.5,
    function: str = 'exponential',
    seed: Optional[int] = None
) -> list:
    """
    Generate random FRB positions within host galaxies.

    FRB positions are sampled from a distribution with
    width proportional to the galaxy's half-light radius.

    Args:
        galaxy_sample: Selected host galaxies
        scale: Scale factor for exponential half-light radius (smaller = more concentrated), or for the
               outer cutoff for the uniform function (arcseconds)
        function: Function to use for generating galaxy positions (exponential, uniform)
        seed: Random seed

    Returns:
        List of SkyCoord objects for true FRB positions
    """
    if seed is not None:
        np.random.seed(seed)

    n_frbs = len(galaxy_sample)

    # Galaxy center coordinates
    galaxy_coords = SkyCoord(
        ra=galaxy_sample.ra.values,
        dec=galaxy_sample.dec.values,
        unit='deg'
    )

    if function == 'exponential':
        # Gamma(2, scale) gives p(r) ∝ r·exp(-r/scale), matching the
        # PATH per-solid-angle exponential prior
        randn = np.random.gamma(shape=2, scale=scale, size=10 * n_frbs)
        good = randn < 6.
        randn = randn[good][:n_frbs]
    elif function == 'exponential_incorrect':
        randn = np.random.exponential(scale=scale, size=10 * n_frbs)
        good = np.abs(randn) < (6.)
        randn = randn[good][:n_frbs]
    elif function == 'uniform_1d':
        # Uniform when integrated over azimuth
        randn = np.random.uniform(low=0., high=10., size=10 * n_frbs)
        good = np.abs(randn) < scale
        randn = randn[good][:n_frbs]
    elif function == 'uniform_2d':
        # scale * sqrt(U) gives p(r) ∝ r over [0, scale], matching the
        # PATH per-solid-angle uniform prior (constant per pixel over a disk)
        randn = scale * np.sqrt(np.random.uniform(low=0., high=1., size=n_frbs))
    else:
        raise ValueError(f"Invalid offset function: {function} (options: 'exponential', 'exponential_incorrect', 'uniform_1d', uniform_2d'")

    if len(randn) < n_frbs:
        raise RuntimeError(
            f"Offset function '{function}' produced only {len(randn)} of {n_frbs} required "
            f"draws after truncation (scale={scale}). Increase the oversampling factor or "
            f"loosen the truncation."
        )

    # Generate offsets
    galaxy_offsets = randn * galaxy_sample.half_light.values * units.arcsec
    position_angles = np.random.uniform(size=n_frbs, low=0., high=360.)

    print("Generating FRB positions within galaxies...")

    # Offset coordinates from galaxy centers (vectorized)
    frb_coords = galaxy_coords.directional_offset_by(
        position_angles * units.deg, galaxy_offsets
    )

    return frb_coords


def _apply_localization_error(
    true_coords,
    localization: np.ndarray,
    seed: Optional[int] = None
):
    """
    Apply per-FRB localization error to true FRB coordinates.

    Offsets are drawn from a truncated (3-sigma) normal distribution
    along each FRB's own error ellipse axes.

    Args:
        true_coords: True FRB coordinates (SkyCoord array or list of SkyCoord)
        localization (np.ndarray): (n_frbs, 3) array of (a, b, PA) in
            (arcsec, arcsec, deg). Produced by `_resolve_localizations`; every
            row may differ.
        seed: Random seed

    Returns:
        Tuple of (observed SkyCoord array, offset magnitudes array in arcsec)
    """
    if seed is not None:
        np.random.seed(seed)

    coords = true_coords if isinstance(true_coords, SkyCoord) else SkyCoord(true_coords)
    n_frbs = len(coords)

    loc = np.asarray(localization, dtype=float)
    if loc.shape != (n_frbs, 3):
        raise ValueError(
            f"localization array must have shape ({n_frbs}, 3), got {loc.shape}. "
            "Pass the output of _resolve_localizations."
        )
    a_ax, b_ax, pa = loc[:, 0], loc[:, 1], loc[:, 2]

    # Generate offsets along major and minor axes (truncated to 3 sigma)
    randn = np.random.normal(size=10*n_frbs)
    # good = np.abs(randn) < 3.
    # randn = randn[good]
    if len(randn) < 2 * n_frbs:
        raise RuntimeError(
            f"Truncated-normal draw produced only {len(randn)} of {2 * n_frbs} required values"
        )

    a_offsets = randn[:n_frbs] * a_ax * units.arcsec
    b_offsets = randn[n_frbs:2 * n_frbs] * b_ax * units.arcsec

    # Total offset magnitude
    loc_offsets = np.sqrt(a_offsets ** 2 + b_offsets ** 2)

    print("Applying localization error...")

    # Apply offsets along each FRB's own ellipse axes (vectorized)
    obs_coords = coords.directional_offset_by(pa * units.deg, a_offsets)
    obs_coords = obs_coords.directional_offset_by((pa + 90.) * units.deg, b_offsets)

    return obs_coords, loc_offsets.value


def _build_output_dataframe(
    obs_coords,
    true_coords,
    galaxy_sample: pd.DataFrame,
    loc_offsets: np.ndarray,
    localization: np.ndarray,
    frb_original_indices: np.ndarray
) -> pd.DataFrame:
    """
    Build output DataFrame with all FRB/host association information.

    Args:
        obs_coords: Observed FRB coordinates (with localization error)
        true_coords: True FRB coordinates (in galaxy)
        galaxy_sample: Assigned host galaxies
        loc_offsets: Localization offset magnitudes (arcsec)
        localization (np.ndarray): (n_frbs, 3) per-FRB (a, b, PA)
        frb_original_indices: Original indices from input FRB DataFrame

    Returns:
        DataFrame with assignment results
    """
    obs = obs_coords if isinstance(obs_coords, SkyCoord) else SkyCoord(obs_coords)
    true = true_coords if isinstance(true_coords, SkyCoord) else SkyCoord(true_coords)

    # Calculate galaxy offsets
    galaxy_coords = SkyCoord(
        ra=galaxy_sample.ra.values,
        dec=galaxy_sample.dec.values,
        unit='deg'
    )

    gal_offsets = true.separation(galaxy_coords).to(units.arcsec).value

    loc = np.asarray(localization, dtype=float)

    df = pd.DataFrame({
        'ra': obs.ra.deg,
        'dec': obs.dec.deg,
        'true_ra': true.ra.deg,
        'true_dec': true.dec.deg,
        'gal_ID': galaxy_sample.ID.values,
        'gal_off': gal_offsets,
        'mag': galaxy_sample.mag.values,
        'half_light': galaxy_sample.half_light.values,
        'loc_off': loc_offsets,
        'FRB_ID': frb_original_indices,
        'a': loc[:, 0],
        'b': loc[:, 1],
        'PA': loc[:, 2],
    })

    return df


def assign_frbs_to_hosts_from_files(
    frb_file: str,
    galaxy_catalog: pd.DataFrame,
    localization: LocalizationLike,
    outfile: str,
    **kwargs
) -> pd.DataFrame:
    """
    Convenience function to assign FRBs to hosts from file and save results.

    This function mirrors the interface of path_simulations.frbs.assign_chime_frbs_to_hosts().

    Args:
        frb_file (str): Path to CSV file containing FRB data (output of generate_frbs)
        galaxy_catalog (pd.DataFrame): Galaxy catalog DataFrame
        localization: (a, b, PA) triple or ellipse population; see
            `assign_frbs_to_hosts`.
        outfile (str): Path to output CSV file
        **kwargs: Additional arguments passed to assign_frbs_to_hosts()

    Returns:
        pd.DataFrame: Assignment results (also saved to outfile)
    """
    # Load FRBs
    frbs = pd.read_csv(frb_file)

    # Assign to hosts
    df = assign_frbs_to_hosts(frbs, galaxy_catalog, localization, **kwargs)

    # Save to disk
    df.to_csv(outfile, index=False)
    print(f"Wrote: {outfile}")

    return df
