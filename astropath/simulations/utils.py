import pandas
import numpy as np
from astropy.coordinates import SkyCoord, match_coordinates_sky
from scipy.signal import convolve
from scipy.special import ellipe
import time


def build_digest(raw_sim_results:pandas.DataFrame=None, frbs:pandas.DataFrame=None, hosts:pandas.DataFrame=None, combined_catalog:pandas.DataFrame=None, 
                 output_fn:str=None, thresh_cross_match:float=2.):
    """
    Combines together the pandas.DataFrame from each simulation step into a "digest"
    DataFrame that can be easily parsed to make informative plots about the simulation results

    Args:
        raw_sim_results (pandas.DataFrame): The simulation results dataframe
        frbs (pandas.DataFrame): The FRBs dataframe
        hosts (pandas.DataFrame): The hosts dataframe
        combined_catalog (pandas.DataFrame): The combined catalog dataframe
        output_fn (str): The filename to save the digest to
        thresh_cross_match (float): A factor that sets the cross-matching threshold
            that determines whether an association is "correct" or not. Specifically,
            it is a multiple of the galaxy half-light radius.

    Returns:
        pandas.DataFrame: The digest dataframe with the following columns:
            - `ra_loc`: Observed FRB RA (degrees) - includes localization error
            - `dec_loc`: Observed FRB Dec (degrees) - includes localization error
            - `true_ra`: True FRB RA in the galaxy (degrees)
            - `true_dec`: True FRB Dec in the galaxy (degrees)
            - `host_ID`: ID of assigned host galaxy within the possible host catalog (see `assign_host` step)
            - `gal_off`: Offset from galaxy center (arcsec)
            - `mag_host`: Apparent r-band magnitude of the true host galaxy
            - `ang_size_host`: Angular size of the true host galaxy (arcsec)
            - `loc_off`: Localization error offset (arcsec)
            - `FRB_ID`: FRB index from the `generate_frbs` step
            - `a`: Localization semi-major axis (arcsec)
            - `b`: Localization semi-minor axis (arcsec)
            - `PA`: Localization position angle (degrees)
            - `ra_host`: RA of the host galaxy center (degrees)
            - `dec_host`: DEC of the host galaxy center (degrees)
            - `sep_best_host_arcsec`: Separation between the center of the best candidate and true host (arcsec)
            - `sep_host_loc_arcsec`: Separation between the center of the true host and localization (arcsec)
            - `sep_best_loc_arcsec`: Separation between the center of the best candidate and localization (arcsec)
            - `sep_host_loc_norm`: Separation between the center of the true host and localization, normalized by the angular size of the true host
            - `sep_best_loc_norm`: Separation between the center of the best candidate and localization, normalized by the angular size of the best candidate
            - `z_host`: Simulated FRB redshift
            - `dmex_host`: Simulated FRB extragalactic DM
            - `frb_mr`: Simulated FRB host apparent r-band magnitude
            - `frb_Mr`: Simulated FRB host absolute r-band magnitude
            - `ra_cand`: RA of the best candidate galaxy center (degrees)
            - `dec_cand`: DEC of the best candidate galaxy center (degrees)
            - `mag_cand`: Apparent r-band magnitude of the best candidate galaxy
            - `ang_size_cand`: Angular size of the best candidate galaxy (arcsec)
            - `cand_ID`: ID of assigned host galaxy within the galaxy catalog used to run PATH (see `run_path` step)
            - `P_O`: Value of the PATH prior P(Oi) for the best candidate
            - `p_xO`: Value of the PATH likelihood p(x|Oi) for the best candidate
            - `P_Ox`: Value of the PATH posterior P(Oi|x) for the best candidate
            - `P_Ux`: Value of the PATH posterior P(U|x)
            - `correct_association`: A boolean indicating whether the best candidate matches the true host (a "correct" association), based on a spatial cross-match
    """

    print("Get parameters from the simulation results dataframe")
    # (like galaxy ID, ra, dec, angular size, magnitude, separation, PATH parameters, etc)
    # and append them to the hosts dataframe to construct a dataframe useful for plotting
    best_cands_list = []
    true_ras = []
    true_decs = []
    true_mags = []
    true_ang_size = []
    true_z = []
    true_dmeg = []
    true_dmhost = []
    true_dmcosmic = []
    true_mr = []
    true_Mr = []
    for ii in range(len(hosts)):
        host_row = hosts[ii:ii+1]
        cands = raw_sim_results[raw_sim_results['iFRB'] == ii]
        frb = frbs.iloc[[ii]]
        # if ii == 106:
        #     print(cands)
        if len(cands) == 0:
            print(ii)
        best_cand = cands[0:1]
        best_cands_list.append(best_cand)

        orig_true_host = combined_catalog[combined_catalog['ID'] == host_row['gal_ID'].item()]
        true_ras.append(orig_true_host.ra.item())
        true_decs.append(orig_true_host.dec.item())
        true_z.append(frb['z'].values[0])
        true_dmeg.append(frb['DMeg'].values[0])
        true_mr.append(frb['m_r'].values[0])
        true_Mr.append(frb['M_r'].values[0])
    best_cands = pandas.concat(best_cands_list, ignore_index=True)

    print("Rename some columns to make concatenation cleaner")
    best_cands = best_cands.rename(
        columns={
            'ra': 'ra_cand', 
            'dec': 'dec_cand',
            'mag': 'mag_cand',
            'ang_size': 'ang_size_cand',
            'ID': 'cand_ID',
        },
    )
    best_cands.drop(columns=['sep', 'iFRB', 'gal_ID'], inplace=True)
    hosts = hosts.rename(
        columns={
            'ra': 'ra_loc', 
            'dec': 'dec_loc',
            'mag': 'mag_host',
            'half_light': 'ang_size_host',
            'gal_ID': 'host_ID',
        },
    )

    print("Calculate angular separation between true host and best candidate")
    true_host_coord = SkyCoord(ra=true_ras, dec=true_decs, unit='deg')
    best_cand_coord = SkyCoord(ra=best_cands.ra_cand.values, dec=best_cands.dec_cand.values, unit='deg')
    sep_best_host_arcsec = true_host_coord.separation(best_cand_coord).arcsec

    
    loc_coord = SkyCoord(ra=hosts.ra_loc.values, dec=hosts.dec_loc.values, unit='deg')
    print("Calculate offset between localization and true host")
    sep_host_loc_arcsec = true_host_coord.separation(loc_coord).arcsec
    sep_host_loc_norm = sep_host_loc_arcsec / hosts.ang_size_host.values
    print("Calculate offset between localization and best candidate")
    sep_best_loc_arcsec = best_cand_coord.separation(loc_coord).arcsec
    sep_best_loc_norm = sep_best_loc_arcsec / best_cands.ang_size_cand.values
    

    print("Add the RA/Dec of the host *galaxy* from the original catalog")
    hosts['ra_host'] = true_ras
    hosts['dec_host'] = true_decs
    hosts['sep_best_host_arcsec'] = sep_best_host_arcsec
    hosts['sep_host_loc_arcsec'] = sep_host_loc_arcsec
    hosts['sep_host_loc_norm'] = sep_host_loc_norm
    hosts['sep_best_loc_arcsec'] = sep_best_loc_arcsec
    hosts['sep_best_loc_norm'] = sep_best_loc_norm
    hosts['z_host'] = true_z
    hosts['dmex_host'] = true_dmeg
    hosts['frb_mr'] = true_mr
    hosts['frb_Mr'] = true_Mr

    print("Merge dataframes, to create a nice big cross-checked dataframe for plotting purposes")
    df = hosts.merge(best_cands, left_index=True, right_index=True)

    print("Determine correct and incorrect matches")
    max_ang_size = np.maximum.reduce([df.ang_size_cand.values, df.ang_size_host.values], axis=0)
    match_criteria = (df.sep_best_host_arcsec.values < thresh_cross_match * max_ang_size)
    df['correct_association'] = match_criteria

    if output_fn is not None:
        print("Saving to file: {}".format(output_fn))
        df.to_parquet(output_fn)
    
    return df


def calculate_unseen(hosts:pandas.DataFrame, galaxy_catalog:pandas.DataFrame,
                     mag_limit:float=15., thresh_cross_match:float=2.):
    """
    Determines whether a set of FRB host galaxies (hosts) is "unseen" (not detected
    as a source) in the given galaxy catalog of limited magnitude depths

    Args:
        hosts (pandas.DataFrame): Simulated FRB hosts, must contain columns:
            `mag_host`, `ra_host`, `dec_host`, `ang_size_host`
        galaxy_catalog (pandas.DataFrame): Galaxy catalog of limited depth,
            must contain columns: `ra`, `dec`, `ang_size`
        mag_limit (float): Lower limit to filter hosts catalog. This exists because we
            make cuts in our galaxy_catalog samples mr < 14. Thus, below this threshold,
            a simulated host will artificially not have a matching source in the galaxy_catalog.
        thresh_cross_match (float): A factor that sets the cross-matching threshold
            that determines whether an FRB host is "unseen" or not. Specifically,
            it is a multiple of the galaxy half-light radius.

    Returns:
        pandas.DataFrame: Hosts catalog, but with an `unseen` column indicating
            whether the host is visible in the galaxy catalog
    """
    if mag_limit is not None:
        hosts = hosts[hosts['mag_host'] > mag_limit]
    
    host_coords = SkyCoord(ra=hosts.ra_host.values, dec=hosts.dec_host.values, unit='deg')
    galaxy_catalog_coords = SkyCoord(ra=galaxy_catalog.ra.values, dec=galaxy_catalog.dec.values, unit='deg')

    # Match hosts to galaxy catalog sources, finding the closest source at any distance
    print("Starting cross-matching (will take a minute)...")
    start_write = time.time()
    idx, d2d, _ = match_coordinates_sky(host_coords, galaxy_catalog_coords, nthneighbor=1)
    end_write = time.time()
    print("Cross-matching finished. Time elapsed: {0:.3f} minutes".format((end_write-start_write)/60.))

    # Filter by distance, finding those matching within the angular size of the galaxy
    max_ang_size = thresh_cross_match * np.maximum.reduce([hosts['ang_size_host'], galaxy_catalog.iloc[idx]['ang_size']], axis=0)
    dup = (d2d.arcsec <  max_ang_size)
    # galaxy_catalog_keep = galaxy_catalog.iloc[idx[dup]]

    # Add "unseen" status to hosts DataFrame
    unseen = np.full(len(hosts), True)
    unseen[dup] = False
    hosts['unseen'] = unseen
    
    return hosts


def azimuthal_integrated_prior(u, theta_prior):
    """
    1D radial PDF p(u) where u = theta/phi, obtained by integrating 
    pw_Oi over azimuth:
        p(u) = 2*pi*theta * pw_Oi(theta) * phi   [change of variables theta->u]

    For 'exp':     p(u) = u*exp(-u/s) / (s^2 * (1-(1+max)*exp(-max)))
    For 'uniform': p(u) = 2u / max^2  for u in [0, max]

    Parameters
    ----------
    u : np.ndarray
        Normalized galactocentric offset theta/phi
    theta_prior : dict
        Same format as pw_Oi: keys 'PDF', 'scale', 'max'

    Returns
    -------
    np.ndarray : normalized 1D PDF values on grid u
    """
    u = np.asarray(u, dtype=float)
    p = np.zeros_like(u)

    if theta_prior['PDF'] == 'exp':
        s     = theta_prior.get('scale', 1.0)
        max_v = theta_prior['max']          # cutoff in units of phi_eff = phi*s
        max_u = max_v * s                   # cutoff in units of phi
        ok    = u <= max_u
        # Normalization: integral of u*exp(-u/s) du from 0 to max_u
        # = s^2 * (1 - (1+max_v)*exp(-max_v))
        norm  = s**2 * (1 - (1 + max_v) * np.exp(-max_v))
        p[ok] = u[ok] * np.exp(-u[ok] / s) / norm

    elif theta_prior['PDF'] == 'uniform':
        max_u = theta_prior['max']
        ok    = u <= max_u
        # p(u) = 2u/max^2, integrates to 1 over [0, max]
        p[ok] = 2 * u[ok] / max_u**2

    elif theta_prior['PDF'] == 'core':
        # p_2d = phi/(theta+phi)/norm  =>  p_1d(u) = 2*pi*phi*u / (u*phi+phi) * phi / norm
        # = 2*pi*phi^2 * u / (phi*(u+1)) / norm = 2*pi*phi * u/(u+1) / norm
        # norm_2d = 2*pi*(term_max - term0), so in normalized units:
        # p(u) = u/(u+1) / integral[u/(u+1) du from 0 to max]
        max_u = theta_prior['max']
        ok    = u <= max_u
        norm  = max_u - np.log(1 + max_u)
        p[ok] = u[ok] / (u[ok] + 1) / norm

    return p


def convolve_prior_with_loc(u, p_prior, sigma_over_phi):
    """
    Convolve the 1D radial prior p(u) (u = theta/phi) with a Gaussian
    localization kernel of width sigma_over_phi = sigma_loc / phi.

    Uses a 1D Gaussian as an approximation to the azimuthally-marginalized
    2D localization error.

    Parameters
    ----------
    u : np.ndarray
        Uniformly spaced grid of theta/phi values
    p_prior : np.ndarray
        Prior values on grid u (from azimuthal_integrated_prior)
    sigma_over_phi : float
        Localization 1-sigma in units of the host half-light radius phi

    Returns
    -------
    u_conv : np.ndarray
        Extended u grid for the convolved distribution
    p_conv : np.ndarray
        Convolved, normalized distribution
    """
    du     = u[1] - u[0]
    kernel = np.exp(-u**2 / (2 * sigma_over_phi**2))
    kernel /= np.sum(kernel) * du          # normalize kernel to unit integral

    p_conv = convolve(p_prior, kernel, mode='full')
    p_conv = np.maximum(p_conv, 0.)        # clip numerical negatives

    u_conv = np.arange(len(p_conv)) * du   # extended grid
    norm   = np.sum(p_conv) * du
    if norm > 0:
        p_conv /= norm

    return u_conv, p_conv


def stack_convolved_prior(u, theta_prior, frb_df, phi_col='ang_size_host'):
    """
    Average the localization-convolved prior over a population of FRBs,
    each with its own localization ellipse and host half-light radius.

    Parameters
    ----------
    u : np.ndarray
        Normalized offset grid (theta/phi)
    theta_prior : dict
        Prior parameters for azimuthal_integrated_prior
    frb_df : pd.DataFrame
        Must have columns 'a', 'b' (loc semi-axes in arcsec) and phi_col
    phi_col : str
        Column name for host half-light radius in arcsec

    Returns
    -------
    u_conv : np.ndarray
    p_stacked : np.ndarray
        Mean convolved prior over all FRBs
    """
    p_prior   = azimuthal_integrated_prior(u, theta_prior)
    stacked   = None

    for _, frb in frb_df.iterrows():
        # Effective circular sigma from ellipse via mean radius of ellipse
        sig_a = max(frb.a, frb.b)
        sig_b = min(frb.a, frb.b)
        k     = np.sqrt(1 - (sig_b / sig_a)**2)
        sigma_arcsec    = (2. / np.pi) * sig_a * ellipe(k) * np.sqrt(2)
        sigma_over_phi  = sigma_arcsec / frb[phi_col]

        u_conv, p_conv = convolve_prior_with_loc(u, p_prior, sigma_over_phi)

        if stacked is None:
            stacked = p_conv.copy()
        else:
            # pad shorter array to match length
            if len(p_conv) < len(stacked):
                p_conv  = np.pad(p_conv,  (0, len(stacked) - len(p_conv)))
            elif len(stacked) < len(p_conv):
                stacked = np.pad(stacked, (0, len(p_conv) - len(stacked)))
            stacked += p_conv

    stacked /= len(frb_df)
    du       = u[1] - u[0]
    stacked /= np.sum(stacked) * du    # renormalize after stacking
    return u_conv, stacked
