import os
import time
import pandas
from astropy.coordinates import SkyCoord

from path_simulations import run
import matplotlib.pyplot as plt
from chime_ffff_pz.path import priors

from IPython import embed

seed = 42
NFRB = 100
output_dir = '/home/xavier/Projects/FRB/data/SIMULATIONS_FINAL/'
survey_str = 'DECaL'

def standard_path_run(catalog_file, hosts_file, output_file='./sim_results_DECaLs.parquet', ncpu=4, return_df=True, save_df=True, multi=True):
    # Get standard DECaLs PATH priors
    prior_dict = priors.load('DECaLS_chime')
    prior_dict['PU'] = 0.15

    final_tbl = run.full(hosts_file, catalog_file, prior_dict, debug=False, ncpu=ncpu, multi=multi)
    
    if save_df:
        final_tbl.to_parquet(output_file)
    
    if return_df:
        return final_tbl

def orig_run_path():

    print("Loading files")
    # Get DECaLs catalog filename to run PATH on
    catalog_fn = os.path.join(output_dir, 'catalog_dudxmmlss_hecate_DECaL_{}_{}.parquet'.format(int(seed), int(NFRB)))
    # Get hosts filename
    hosts_fn = os.path.join(output_dir, 'generated_hosts_DECaL_DECaLhost_hecatecut_{}_{}.parquet'.format(int(seed), int(NFRB)))
    # Set filename of output simulation results
    output_fn = os.path.join(output_dir, 'sim_results_DECaL_DECaLhost_hecatecut_{}_{}.parquet'.format(int(seed), int(NFRB)))

    start = time.time()
    sim_results = standard_path_run(
        catalog_fn, hosts_fn, 
        output_file=output_fn, ncpu=4,
        multi=True,
    )
    end = time.time()
    print("Reading finished. Time elapsed: {0:.3f} seconds".format((end-start)))

def write_digest():

    generated_frbs_fn = os.path.join(output_dir, 'generated_frbs_test_{}_{}.parquet'.format(int(seed), int(NFRB)))
    frbs = pandas.read_parquet(generated_frbs_fn)

    hosts_fn = os.path.join(output_dir, 'generated_hosts_DECaL_DECaLhost_hecatecut_{}_{}.parquet'.format(int(seed), int(NFRB)))
    hosts = pandas.read_parquet(hosts_fn).reset_index()

    #combined_catalog_fn = os.path.join(os.getenv('CHIME_SANDBOX'), '../../', 'catalogs', 'combined_HSC_DECaLs_HECATE_galaxies_hecatecut.parquet')
    #combined_catalog = pandas.read_parquet(combined_catalog_fn)
    combined_file = os.path.join(os.getenv('FRB_APATH'), 'combined_HSC_DECaLs_HECATE_galaxies_hecatecut.parquet')
    combined_catalog = pandas.read_parquet(combined_file)


    # Set the survey string for saving files and making plots
    sim_results_fn = os.path.join(output_dir, 'sim_results_DECaL_DECaLhost_hecatecut_{}_{}.parquet'.format(int(seed), int(NFRB)))
    raw_sim_results = pandas.read_parquet(sim_results_fn)

    print("Get parameters from the simulation results dataframe")
    # (like galaxy ID, ra, dec, angular size, magnitude, separation, PATH parameters, etc)
    # and append them to the hosts dataframe to construct a dataframe useful for plotting
    best_cands_list = []
    true_ras = []
    true_decs = []
    true_mags = []
    true_ang_size = []
    true_z = []
    true_dmex = []
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
        true_dmex.append(frb['DMex'].values[0])
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
    hosts['dmex_host'] = true_dmex
    hosts['frb_mr'] = true_mr
    hosts['frb_Mr'] = true_Mr

    print("Merge dataframes, to create a nice big cross-checked dataframe for plotting purposes")
    df = hosts.merge(best_cands, left_index=True, right_index=True)

    plotting_fn = os.path.join(output_dir, f'plot_data_{survey_str}_centroidfix_zinfo.parquet')
    print("Saving to file: {}".format(plotting_fn))
    df.to_parquet(plotting_fn)

def orig_final_plot():
    plotting_fn = os.path.join(output_dir, f'plot_data_{survey_str}_centroidfix_zinfo.parquet')
    df = pandas.read_parquet(plotting_fn)
    #print(df.head())

    fig = plt.figure(figsize=(8,5))

    match_criteria = df['host_ID'] == df['cand_ID']

    correct_associations = df[match_criteria]
    incorrect_associations = df[~match_criteria]

    # (left, bottom, width, height)
    correct = (0, 0., 0.5, 1.)
    incorrect = (0.5, 0., 0.5, 1.)


    cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    ax = plt.axes(correct)
    plt.scatter(correct_associations['mag_cand'], correct_associations['P_Ox'], color=cycle[0], marker='o', s=30, alpha=0.5, label='Correct Association')
    plt.xlabel(r'$m_r$: Best Candidate')
    plt.ylabel(r'PATH $P(O|x)$: Best Candidate')
    plt.title('True Associations: {}'.format(len(correct_associations)))
    plt.grid(zorder=20)
    plt.ylim([0,1])
    plt.xlim([10,29])
    #     plt.axvline(x=18, color='black')
    plt.axhline(y=0.9, color='black')

    ax = plt.axes(incorrect)
    plt.scatter(incorrect_associations['mag_cand'], incorrect_associations['P_Ox'], color=cycle[1], marker='x', s=30, alpha=0.5, label='Incorrect Association')
    plt.xlabel(r'$m_r$: Best Candidate')
    ax.set_yticklabels([])
    plt.title('False Associations: {}'.format(len(incorrect_associations)))
    plt.grid(zorder=20)
    plt.ylim([0,1])
    plt.xlim([10,29])
    #     plt.axvline(x=18, color='black')
    plt.axhline(y=0.9, color='black')

    #plt.show()
    plt.savefig(os.path.join(output_dir, f'plot_mr_POx_{survey_str}_centroidfix_zinfo.png'), dpi=200, format="png", bbox_inches="tight")
    embed(header='151 of orig_final_plot')

# Command line
if __name__ == "__main__":
    #orig_run_path()
    #write_digest()
    orig_final_plot()