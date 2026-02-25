import os
import time
import numpy as np
import pandas
from astropy.coordinates import SkyCoord

from path_simulations import run as orig_run
import matplotlib.pyplot as plt

from astropath.simulations import run_path
from astropath.simulations import utils as sim_utils

from IPython import embed

# Parameters
seed = 42
NFRB = 100
output_dir = '/home/xavier/Projects/FRB/data/SIMULATIONS_FINAL/'
survey_str = 'DECaL'

# Standard DECaL Priors
DECaLS_chime_priors = {}
DECaLS_chime_priors['version'] = '1.0'
DECaLS_chime_priors['PU'] = 0.15
DECaLS_chime_priors['survey'] = 'DECaL'  # No S in the FRB repo..
DECaLS_chime_priors['scale'] = 0.5
prior_dict = DECaLS_chime_priors

def standard_path_run(catalog_file, hosts_file, output_file='./sim_results_DECaLs.parquet', ncpu=4, return_df=True, save_df=True, multi=True):

    final_tbl = orig_run.full(hosts_file, catalog_file, prior_dict, debug=False, ncpu=ncpu, multi=multi)
    
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

def build_orig_digest():

    generated_frbs_fn = os.path.join(output_dir, 'generated_frbs_test_{}_{}.parquet'.format(int(seed), int(NFRB)))
    frbs = pandas.read_parquet(generated_frbs_fn)
    frbs['DMeg'] = frbs['DMex']

    hosts_fn = os.path.join(output_dir, 'generated_hosts_DECaL_DECaLhost_hecatecut_{}_{}.parquet'.format(int(seed), int(NFRB)))
    hosts = pandas.read_parquet(hosts_fn).reset_index()

    #combined_catalog_fn = os.path.join(os.getenv('CHIME_SANDBOX'), '../../', 'catalogs', 'combined_HSC_DECaLs_HECATE_galaxies_hecatecut.parquet')
    #combined_catalog = pandas.read_parquet(combined_catalog_fn)
    combined_file = os.path.join(os.getenv('FRB_APATH'), 'combined_HSC_DECaLs_HECATE_galaxies_hecatecut.parquet')
    combined_catalog = pandas.read_parquet(combined_file)

    sim_results_fn = os.path.join(output_dir, 'sim_results_DECaL_DECaLhost_hecatecut_{}_{}.parquet'.format(int(seed), int(NFRB)))
    raw_sim_results = pandas.read_parquet(sim_results_fn)

    plotting_fn = os.path.join(output_dir, f'plot_data_{survey_str}_centroidfix_zinfo.parquet')

    digest = sim_utils.build_digest(raw_sim_results, frbs, hosts, combined_catalog, plotting_fn)

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

def test_astropath_path():
    # Get DECaLs catalog filename to run PATH on
    catalog_fn = os.path.join(output_dir, 'catalog_dudxmmlss_hecate_DECaL_{}_{}.parquet'.format(int(seed), int(NFRB)))
    catalog = pandas.read_parquet(catalog_fn)
    # Get hosts filename
    hosts_fn = os.path.join(output_dir, 'generated_hosts_DECaL_DECaLhost_hecatecut_{}_{}.parquet'.format(int(seed), int(NFRB)))
    hosts = pandas.read_parquet(hosts_fn).reset_index()
    # Set filename of output simulation results
    #output_fn = os.path.join(output_dir, 'sim_results_DECaL_DECaLhost_hecatecut_{}_{}.parquet'.format(int(seed), int(NFRB)))

    final_tbl = run_path.full(hosts, catalog, prior_dict, 
        debug=False, ncpu=4, multi=True)

    # Digest
    generated_frbs_fn = os.path.join(output_dir, 'generated_frbs_test_{}_{}.parquet'.format(int(seed), int(NFRB)))
    frbs = pandas.read_parquet(generated_frbs_fn)
    frbs['DMeg'] = frbs['DMex']
    combined_file = os.path.join(os.getenv('FRB_APATH'), 'combined_HSC_DECaLs_HECATE_galaxies_hecatecut.parquet')
    combined_catalog = pandas.read_parquet(combined_file)

    digest = sim_utils.build_digest(final_tbl, frbs, hosts, combined_catalog)

    # Orig
    plotting_fn = os.path.join(output_dir, f'plot_data_{survey_str}_centroidfix_zinfo.parquet')
    orig_digest = pandas.read_parquet(plotting_fn)

    # Test
    assert np.allclose(digest['host_ID'].values, orig_digest['host_ID'].values)
    assert np.allclose(digest['cand_ID'].values, orig_digest['cand_ID'].values)
    assert np.allclose(digest['ra_host'].values, orig_digest['ra_host'].values)
    assert np.allclose(digest['dec_host'].values, orig_digest['dec_host'].values)
    assert np.allclose(digest['mag_true_host'].values, orig_digest['mag_true_host'].values)
    assert np.allclose(digest['ang_size_host'].values, orig_digest['ang_size_host'].values)
    assert np.allclose(digest['sep_best_host'].values, orig_digest['sep_best_host'].values)
    assert np.allclose(digest['z_host'].values, orig_digest['z_host'].values)
    assert np.allclose(digest['dmex_host'].values, orig_digest['dmex_host'].values)
    assert np.allclose(digest['frb_mr'].values, orig_digest['frb_mr'].values)
    assert np.allclose(digest['frb_Mr'].values, orig_digest['frb_Mr'].values)

    print("Tests passed!!")

# Command line
if __name__ == "__main__":
    #orig_run_path()
    #build_orig_digest()
    #orig_final_plot()

    test_astropath_path()