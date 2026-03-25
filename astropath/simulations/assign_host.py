"""
Module for assigning simulated FRBs to host galaxies by apparent magnitude (m_r).

This module implements magnitude-based host assignment that matches FRBs to
galaxies from a catalog by their apparent r-band magnitudes, ensuring realistic
associations where brighter FRBs tend to be assigned to brighter galaxies.
"""

import os
import numpy as np
import random
import pandas as pd
from pathlib import Path
from typing import Tuple, Optional

from astropy import units
from astropy.coordinates import SkyCoord, match_coordinates_sky

from IPython import embed

def load_galaxy_catalog():
    # Try to load real catalog
    frb_apath = os.environ.get('FRB_APATH')
    
    if frb_apath is not None:
        catalog_path = Path(frb_apath) / 'combined_HSC_DECaLs_HECATE_galaxies_hecatecut.parquet'
        
        if catalog_path.exists():
            print(f"Loading real galaxy catalog from:")
            print(f"  {catalog_path}")
            galaxies = pd.read_parquet(catalog_path)
            
            return galaxies
        else:
            raise ValueError(f"Catalog not found at {catalog_path}")
    else:
        raise ValueError("FRB_APATH not set")



def assign_frbs_to_hosts(
    frb_df: pd.DataFrame,
    galaxy_catalog: pd.DataFrame,
    localization: Tuple[float, float, float],
    mag_range: Tuple[float, float] = None, #(17., 28.),
    scale: float = 0.5,
    trim_catalog: units.Quantity = 1 * units.arcmin,
    seed: Optional[int] = None,
    debug: bool = False
) -> pd.DataFrame:
    """
    Assign FRBs to host galaxies based on apparent magnitude matching.

    This function assigns each FRB to a host galaxy by matching their apparent
    magnitudes (m_r). The matching algorithm uses a "fake coordinate" approach
    where magnitudes are encoded as declinations, allowing sky coordinate
    matching to effectively match by brightness.

    Each FRB is randomly placed within the host galaxy according to an
    exponential offset distribution, then the observed coordinates are offset
    according to the localization error ellipse.

    Args:
        frb_df (pd.DataFrame): FRB catalog with columns:
            - 'm_r': Apparent r-band magnitude of the host
            Additional columns are preserved in output
        galaxy_catalog (pd.DataFrame): Galaxy catalog with columns:
            - 'ra': Right ascension (degrees)
            - 'dec': Declination (degrees)
            - 'mag_best': Apparent r-band magnitude
            - 'half_light': Half-light radius (arcsec)
            - 'ID': Unique galaxy identifier
        localization (tuple): Error ellipse parameters (a, b, PA) where:
            - a: Semi-major axis (arcsec)
            - b: Semi-minor axis (arcsec)
            - PA: Position angle (degrees, East of North)
        mag_range (tuple, optional): (min, max) magnitude range for FRB selection.
            FRBs outside this range are filtered out. Default: (17., 28.)
        scale (float, optional): Scale factor for galaxy half-light radius when
            placing FRBs. Smaller values concentrate FRBs closer to galaxy centers.
            Default: 0.5
        trim_catalog (units.Quantity, optional): Buffer to trim from catalog edges
            to ensure FRBs stay within analysis region. Default: 1 arcmin
        seed (int, optional): Random seed for reproducibility
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
            - 'a': Localization semi-major axis (arcsec)
            - 'b': Localization semi-minor axis (arcsec)
            - 'PA': Localization position angle (degrees)

    Raises:
        ValueError: If required columns are missing from input DataFrames

    Example:
        >>> frbs = generate_frbs(1000, 'CHIME')
        >>> # Load galaxy catalog
        >>> galaxies = pd.read_csv('galaxy_catalog.csv')
        >>> # Assign with 0.5" x 0.3" error ellipse at PA=45deg
        >>> assignments = assign_frbs_to_hosts(
        ...     frbs, galaxies, localization=(0.5, 0.3, 45)
        ... )
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

    # Trim catalog edges to maintain analysis region
    galaxy_cut = _trim_catalog(galaxy_catalog, trim_catalog)

    if len(galaxy_cut) == 0:
        raise ValueError("No galaxies remain after trimming catalog edges")

    # Match FRBs to galaxies by magnitude
    galaxy_indices = _match_by_magnitude(cut_frbs, galaxy_cut, debug=debug)
    galaxy_sample = galaxy_cut.loc[galaxy_indices]

    # Generate FRB positions within galaxies
    true_coords = _generate_galaxy_positions(
        galaxy_sample, scale=scale, #seed=seed
    )

    # Apply localization error
    obs_coords, loc_offsets = _apply_localization_error(
        true_coords, localization, #seed=seed
    )

    # Build output DataFrame
    df_out = _build_output_dataframe(
        obs_coords, true_coords, galaxy_sample,
        loc_offsets, localization, cut_frbs.index.values
    )

    return df_out


def _validate_frb_columns(frb_df: pd.DataFrame):
    """Validate that FRB DataFrame has required columns."""
    required_cols = ['m_r']
    missing = [col for col in required_cols if col not in frb_df.columns]
    if missing:
        raise ValueError(f"FRB DataFrame missing required columns: {missing}")


def _validate_galaxy_columns(galaxy_df: pd.DataFrame):
    """Validate that galaxy DataFrame has required columns."""
    required_cols = ['ra', 'dec', 'mag_best', 'half_light', 'ID']
    missing = [col for col in required_cols if col not in galaxy_df.columns]
    if missing:
        raise ValueError(f"Galaxy catalog missing required columns: {missing}")


def _trim_catalog(galaxy_df: pd.DataFrame, trim: units.Quantity) -> pd.DataFrame:
    """
    Trim edges off catalog to maintain PATH analysis region.

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

    This implements the clever magnitude-matching algorithm from
    frb.frb_surveys.mock.frbs_in_hosts(). It creates "fake" sky coordinates
    where the declination encodes the magnitude, then uses astropy's
    match_coordinates_sky to match by brightness.

    The algorithm iteratively matches FRBs to galaxies, ensuring each galaxy
    is used only once. Brighter FRBs are preferentially matched to brighter
    galaxies.

    Args:
        frb_df: FRB catalog with 'm_r' column
        galaxy_df: Galaxy catalog with 'mag_best' column
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
        dec=galaxy_df['mag_best'].values,
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
    function:str='exponential',
    seed: Optional[int] = None
) -> list:
    """
    Generate random FRB positions within host galaxies.

    FRB positions are sampled from a truncated normal distribution with
    width proportional to the galaxy's half-light radius.

    Args:
        galaxy_sample: Selected host galaxies
        scale: Scale factor for half-light radius (smaller = more concentrated)
        function: Function to use for generating galaxy positions (exponential, uniform, truncated normal)
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

    # Generate offsets from galaxy centers
    # Use truncated normal distribution (within 6 sigma)
    #theta_max = galaxy_sample.half_light.values / scale
    #randn = np.random.normal(size=10 * n_frbs)
    #good = np.abs(randn) < (6. * scale)
    #randn = randn[good][:n_frbs]

    if function == 'exponential':
        randn = np.random.exponential(scale=scale, size=10 * n_frbs)
        good = np.abs(randn) < (6.)
        randn = randn[good][:n_frbs]
    elif function == 'uniform':
        randn = np.random.uniform(low=0., high=10., size=10*n_frbs)
        good = np.abs(randn) < (6.)
        randn = randn[good][:n_frbs]
    #elif function == 'truncated normal':
    #    randn = np.random.normal(size=10 * n_frbs)
    #    good = np.abs(randn) < (6.)
    #    randn = randn[good][:n_frbs]
    else:
        raise ValueError(f"Invalid function: {function}")

    # Generate offsets
    galaxy_offsets = randn * galaxy_sample.half_light.values * units.arcsec
    position_angles = np.random.uniform(size=n_frbs, low=0., high=360.)
    #embed(header='assign_host 382')

    print("Generating FRB positions within galaxies...")

    # Offset coordinates from galaxy centers
    frb_coords = [
        coord.directional_offset_by(position_angles[kk] * units.deg, galaxy_offsets[kk])
        for kk, coord in enumerate(galaxy_coords)
    ]

    return frb_coords


def _apply_localization_error(
    true_coords: list,
    localization: Tuple[float, float, float],
    seed: Optional[int] = None
) -> Tuple[list, np.ndarray]:
    """
    Apply localization error to true FRB coordinates.

    Offsets are drawn from a truncated (3-sigma) normal distribution
    along the error ellipse axes.

    Args:
        true_coords: List of true FRB SkyCoord objects
        localization: (a, b, PA) error ellipse parameters (arcsec, arcsec, deg)
        seed: Random seed

    Returns:
        Tuple of (observed coordinates list, offset magnitudes array)
    """
    if seed is not None:
        np.random.seed(seed)

    n_frbs = len(true_coords)
    a, b, pa = localization

    # Generate offsets along major and minor axes (truncated to 3 sigma)
    randn = np.random.normal(size=10 * n_frbs)
    good = np.abs(randn) < 3.
    randn = randn[good]

    a_offsets = randn[:n_frbs] * a * units.arcsec
    b_offsets = randn[n_frbs:2*n_frbs] * b * units.arcsec

    # Total offset magnitude
    loc_offsets = np.sqrt(a_offsets**2 + b_offsets**2)

    print("Applying localization error...")

    # Apply offsets along ellipse axes
    obs_coords = [
        coord.directional_offset_by(pa * units.deg, a_offsets[kk])
        for kk, coord in enumerate(true_coords)
    ]
    obs_coords = [
        coord.directional_offset_by((pa + 90) * units.deg, b_offsets[kk])
        for kk, coord in enumerate(obs_coords)
    ]

    return obs_coords, loc_offsets.value


def _build_output_dataframe(
    obs_coords: list,
    true_coords: list,
    galaxy_sample: pd.DataFrame,
    loc_offsets: np.ndarray,
    localization: Tuple[float, float, float],
    frb_original_indices: np.ndarray
) -> pd.DataFrame:
    """
    Build output DataFrame with all FRB/host association information.

    Args:
        obs_coords: Observed FRB coordinates (with localization error)
        true_coords: True FRB coordinates (in galaxy)
        galaxy_sample: Assigned host galaxies
        loc_offsets: Localization offset magnitudes (arcsec)
        localization: Error ellipse parameters
        frb_original_indices: Original indices from input FRB DataFrame

    Returns:
        DataFrame with assignment results
    """
    # Calculate galaxy offsets
    galaxy_coords = SkyCoord(
        ra=galaxy_sample.ra.values,
        dec=galaxy_sample.dec.values,
        unit='deg'
    )

    gal_offsets = [
        coord.separation(galaxy_coords[kk]).to(units.arcsec).value
        for kk, coord in enumerate(true_coords)
    ]

    df = pd.DataFrame({
        'ra': [coord.ra.deg for coord in obs_coords],
        'dec': [coord.dec.deg for coord in obs_coords],
        'true_ra': [coord.ra.deg for coord in true_coords],
        'true_dec': [coord.dec.deg for coord in true_coords],
        'gal_ID': galaxy_sample.ID.values,
        'gal_off': gal_offsets,
        'mag': galaxy_sample.mag_best.values,
        'half_light': galaxy_sample.half_light.values,
        'loc_off': loc_offsets,
        'FRB_ID': frb_original_indices,
        'a': localization[0],
        'b': localization[1],
        'PA': localization[2]
    })

    return df


def assign_frbs_to_hosts_from_files(
    frb_file: str,
    galaxy_catalog: pd.DataFrame,
    localization: Tuple[float, float, float],
    outfile: str,
    **kwargs
) -> pd.DataFrame:
    """
    Convenience function to assign FRBs to hosts from file and save results.

    This function mirrors the interface of path_simulations.frbs.assign_chime_frbs_to_hosts().

    Args:
        frb_file (str): Path to CSV file containing FRB data (output of generate_frbs)
        galaxy_catalog (pd.DataFrame): Galaxy catalog DataFrame
        localization (tuple): (a, b, PA) error ellipse parameters (arcsec, arcsec, deg)
        outfile (str): Path to output CSV file
        **kwargs: Additional arguments passed to assign_frbs_to_hosts()

    Returns:
        pd.DataFrame: Assignment results (also saved to outfile)

    Example:
        >>> df = assign_frbs_to_hosts_from_files(
        ...     'frbs_chime.csv',
        ...     galaxy_catalog,
        ...     localization=(0.5, 0.3, 45),
        ...     outfile='frb_host_assignments.csv'
        ... )
    """
    # Load FRBs
    frbs = pd.read_csv(frb_file)

    # Assign to hosts
    df = assign_frbs_to_hosts(frbs, galaxy_catalog, localization, **kwargs)

    # Save to disk
    df.to_csv(outfile, index=False)
    print(f"Wrote: {outfile}")

    return df
