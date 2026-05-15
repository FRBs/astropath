import pandas
import numpy as np
from astropy.coordinates import SkyCoord


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
    match_criteria = (df.sep_best_host_arcsec.values < thresh_cross_match*max_ang_size)
    df['correct_association'] = match_criteria

    if output_fn is not None:
        print("Saving to file: {}".format(output_fn))
        df.to_parquet(output_fn)
    
    return df
