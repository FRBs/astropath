import pandas
from astropy.coordinates import SkyCoord


def build_digest(raw_sim_results:pandas.DataFrame=None, frbs:pandas.DataFrame=None, hosts:pandas.DataFrame=None, combined_catalog:pandas.DataFrame=None, 
    output_fn:str=None):
    """
    Construct a digest of the simulation results.

    Args:
        raw_sim_results (pandas.DataFrame): The simulation results dataframe
        frbs (pandas.DataFrame): The FRBs dataframe
        hosts (pandas.DataFrame): The hosts dataframe
        combined_catalog (pandas.DataFrame): The combined catalog dataframe
        output_fn (str): The filename to save the digest to

    Returns:
        pandas.DataFrame: The digest dataframe
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
        true_mags.append(orig_true_host.mag_best.item())
        true_ang_size.append(orig_true_host.ang_size.item())
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
            'sep': 'sep_cand',
            'ID': 'cand_ID',
            'gal_ID': 'gal_ind',
        },
    )
    hosts = hosts.rename(
        columns={
            'ra': 'ra_loc', 
            'dec': 'dec_loc',
            'mag': 'mag_host',
            'half_light': 'half_light_host',
            'gal_ID': 'host_ID',
        },
    )

    print("Calculate angular separation between true host and best candidate")
    true_host_coord = SkyCoord(ra=true_ras, dec=true_decs, unit='deg')
    best_cand_coord = SkyCoord(ra=best_cands.ra_cand.values, dec=best_cands.dec_cand.values, unit='deg')
    sep = true_host_coord.separation(best_cand_coord).arcsec

    print("Add the RA/Dec of the host *galaxy* from the original catalog")
    hosts['ra_host'] = true_ras
    hosts['dec_host'] = true_decs
    hosts['mag_true_host'] = true_mags
    hosts['ang_size_host'] = true_ang_size
    hosts['sep_best_host'] = sep
    hosts['z_host'] = true_z
    hosts['dmex_host'] = true_dmeg
    hosts['frb_mr'] = true_mr
    hosts['frb_Mr'] = true_Mr

    print("Merge dataframes, to create a nice big cross-checked dataframe for plotting purposes")
    df = hosts.merge(best_cands, left_index=True, right_index=True)

    if output_fn is not None:
        print("Saving to file: {}".format(output_fn))
        df.to_parquet(output_fn)
    
    return df