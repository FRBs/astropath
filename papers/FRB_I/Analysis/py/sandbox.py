""" Routines Analyze the sandbox(es) """
import sys, os
import numpy as np
import glob

import multiprocessing

from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages

import seaborn as sns

import pandas

from astropy.coordinates import SkyCoord
from astropy.coordinates import search_around_sky
from astropy import units

from frb.associate import frbassociate
from frb import frb


# Local
sys.path.append(os.path.abspath("../Analysis/py"))
import associate_defs

from IPython import embed


def cosmos_defs():
    cdefs = dict(plate_scale=0.05,
                 filter='ACS_i',
                 )
    return cdefs


def main_calc(idx, frb_coord, frb_ee, close_galaxies, sb_defs, prior,
              max_radius):
    """

    Args:
        idx:
        frb_coord:
        frb_ee:
        sb_galaxies:
        galaxy_coords:
        sb_defs:
        prior:
        max_radius:

    Returns:
        tuple:  pandas.DataFrame(), int

    """
    # Empty?
    if len(close_galaxies) == 0:
        return None, idx

    # Dummy FRB (to avoid re-instantiation)
    FRB = frb.FRB('sandbox', frb_coord, 1. * units.pc / units.cm ** 3)
    FRB.eellipse = frb_ee.copy()

    #
    frbA = frbassociate.FRBAssociate(FRB)
    frbA.filter = sb_defs['filter']
    frbA.max_radius = max_radius

    # Grab galaxies
    #sep = frb_coord.separation(galaxy_coords)
    #galaxy_cut = (sep < 15 * units.arcsec)
    #galaxy_cut = (sep < 15 * units.arcsec) & (sb_galaxies.class_star < 0.6) & (
    #        sb_galaxies.flux_radius > 1.)

    # Cut
    #close_galaxies = sb_galaxies[galaxy_cut][['mag_best', 'kron_radius', 'a_image',
    #                                          'b_image', 'ra', 'dec', 'class_star']]

    # Update FRB object
    FRB.coord = frb_coord
    frbA.candidates = close_galaxies
    #
    #frbA.candidates['half_light'] = close_galaxies.a_image * sb_defs['plate_scale']
    #frbA.candidates['separation'] = sep[galaxy_cut].to('arcsec').value
    #frbA.candidates[frbA.filter] = close_galaxies.mag_best
    #frbA.candidates['coords'] = galaxy_coords[galaxy_cut]

    # Pchance
    frbA.calc_pchance()

    # Priors
    frbA.calc_priors(prior['U'], method=prior['O'])
    prior['theta']['ang_size'] = frbA.candidates['half_light'].values
    frbA.set_theta_prior(prior['theta'])

    # Posteriors
    frbA.calc_POx()

    # Save
    sv_tbl = frbA.candidates.sort_values('P_Ox', ascending=False).drop(columns=['coords'])
    return sv_tbl, idx, frbA


def run_bayes(frb_tbl_file, galaxy_cat_file, prior, outfile,
              frb_ee,
              galaxy_coords = None,
              max_radius=20.,
              mag_lim=None,
              ncpu=15, multi=True):
    """
    Run the Bayes analysis

    Args:
        frb_tbl_file:
        galaxy_cat_file:
        prior:
        outfile:
        frb_ee:
        galaxy_coords:
        max_radius:

    Returns:

    """
    # Defs
    sb_defs = cosmos_defs()

    # Load FRBs
    print("Loading...")
    sb_frbs = pandas.read_csv(frb_tbl_file, index_col=0)
    nFRB = len(sb_frbs)

    # frb_ee
    if isinstance(frb_ee, dict):
        frb_ee_list = [frb_ee]*nFRB
    elif isinstance(frb_ee, str):
        frb_ee_list = []
        if frb_ee == 'loc_sig':
            for kk in range(nFRB):
                frb_ee_list.append({'a': sb_frbs.iloc[kk].loc_sig,
                                    'b': sb_frbs.iloc[kk].loc_sig,
                                    'theta': 0.})
    else:
        raise IOError("Bad frb_ee")

    # Load Galaxies
    sb_galaxies = pandas.read_feather(galaxy_cat_file)
    # Cut
    sb_galaxies = sb_galaxies[np.isfinite(sb_galaxies.a_image) &
                              np.isfinite(sb_galaxies.mag_best)]

    # Extra bits
    sb_galaxies['half_light'] = sb_galaxies.a_image * sb_defs['plate_scale']
    sb_galaxies[sb_defs['filter']] = sb_galaxies.mag_best

    # Cut on mag?
    if mag_lim is not None:
        keep = sb_galaxies.mag_best < mag_lim
        sb_galaxies = sb_galaxies[keep]

    # Coordinates (this is slow)
    print("Building coords...")
    frb_coords = SkyCoord(ra=sb_frbs.frb_ra, dec=sb_frbs.frb_dec, unit='deg')
    if galaxy_coords is None:
        galaxy_coords = SkyCoord(ra=sb_galaxies.ra, dec=sb_galaxies.dec, unit='deg')

    # Find the candidates
    print("Finding candidates...")
    idx1, idx2, sep2d, _ = search_around_sky(galaxy_coords, frb_coords, 15*units.arcsec)

    # Slice the tables
    print("Slicing...")
    list_candidates = []
    #nFRB = 2000
    for kk in range(nFRB):
        if (kk % 1000) == 0:
            print('kk: ', kk)
        in_idx2 = np.where(idx2 == kk)[0]
        gd_gal = idx1[in_idx2]
        close_galaxies = sb_galaxies.iloc[gd_gal][[#'mag_best', #'a_image', #'kron_radius',
                                                  #'b_image', #'ra', 'dec', 'class_star',
            'half_light', sb_defs['filter'],
            'ra', 'dec']]
        # Extras
        close_galaxies['separation'] = sep2d[gd_gal].to('arcsec').value
        close_galaxies['coords'] = galaxy_coords[gd_gal]
        #
        list_candidates.append(close_galaxies)

    # Loop me!
    print('Starting to loop!')
    # Build
    if multi:
        pool = multiprocessing.Pool(processes=ncpu)
        results = [pool.apply_async(main_calc,
            args=(idx_FRB, frb_coords[idx_FRB], frb_ee_list[idx_FRB],
                    list_candidates[idx_FRB], sb_defs, prior, max_radius))
            for idx_FRB in range(nFRB)]
        output = [p.get() for p in results]
        idx = [item[1] for item in output]
        all_tbls = np.array([item[0] for item in output], dtype=object)
        all_tbls = all_tbls[idx]
        # Expunge the None's
        gd_tbl = np.array([False if item is None else True for item in all_tbls])
        gd_idx = np.arange(all_tbls.size)[gd_tbl]
        all_tbls = all_tbls[gd_tbl]
        for kk in range(all_tbls.size):
            all_tbls[kk]['iFRB'] = gd_idx[kk]

    else:
        embed(header='broken')
        all_tbls = []
        for idx_FRB in range(nFRB):
            sv_tbl = main_calc(frb_coords[idx_FRB], sb_galaxies, galaxy_coords,
                      FRB, frbA, sb_defs, prior)

            sv_tbl['iFRB'] = idx_FRB
            all_tbls.append(sv_tbl)

    # Finish
    final_tbl = pandas.concat(all_tbls)
    final_tbl.to_csv(outfile)
    print("Wrote: {}".format(outfile))


def best_stat_table(results_tbl_file, outfile, lbl='U'):
    """

    Args:
        results_tbl_file:
        outfile:
        lbl:

    Returns:

    """
    # Load
    results = pandas.read_csv(results_tbl_file, index_col=0)

    # Sub tables
    sub_tbls = []
    for uni_iFRB in np.unique(results.iFRB):
        sub_tbls.append(results[results.iFRB == uni_iFRB].sort_values('P_Ox', ascending=False))

    #
    max_POx, best_offset, best_half, best_mag = [], [], [], []
    best_ra, best_dec = [], []
    #
    for sub_tbl in sub_tbls:
        # Max
        max_cand = sub_tbl.iloc[0]
        # Max POx
        max_POx.append(max_cand.P_Ox)
        # Best stuff
        best_offset.append(max_cand.separation)
        best_half.append(max_cand.half_light)
        best_mag.append(max_cand.ACS_i)
        best_ra.append(max_cand.ra)
        best_dec.append(max_cand.dec)
        #
    stats = pandas.DataFrame(dict(max_POx=max_POx, best_offset=best_offset,
                                  best_half=best_half,
                                  best_mag=best_mag,
                                  best_ra=best_ra,
                                  best_dec=best_dec))
    stats['PriorSet'] = lbl
    #
    stats.to_csv(outfile)
    print("Wrote: {}".format(outfile))


def build_pdfs(results_tbl_file, outfile):
    """
    Genrate PDFs of the main quantities

    Args:
        results_tbl_file:
        outfile:

    Returns:

    """

    # Load
    results = pandas.read_csv(results_tbl_file, index_col=0)

    # Do it
    scale = 1000
    model_mags, model_theta = [], []
    # Record
    for cc, cand in results.iterrows():
        N = int(np.round(cand.P_Ox * scale))
        if N == 0:
            continue
        model_mags += [cand.ACS_i] * N
        model_theta += [(cand.separation / cand.half_light)] * N
    # Pandas
    df = pandas.DataFrame(dict(pdf_mags=model_mags, pdf_theta=model_theta))
    # Write
    df.to_csv(outfile)
    print("Wrote: {}".format(outfile))



def do_stats(results_file, stats_outfile):
    """
    Run stats on the results

    Args:
        results_file:
        stats_outfile:

    Returns:

    """
    # Load
    print("Loading results..")
    results = pandas.read_csv(results_file, index_col=0)

    # Cut em up
    print("Cutting..")
    sub_tbls = []
    for uni_iFRB in np.unique(results.iFRB):
        sub_tbls.append(results[results.iFRB == uni_iFRB].sort_values('P_Ox', ascending=False))

    # Stats
    max_POx, best_offset, best_half, best_mag = [], [], [], []
    best_ra, best_dec, best_Pc, best_PO = [], [], [], []
    PU, iFRB = [], []
    print("Stats..")
    #
    for sub_tbl in sub_tbls:
        # Max
        max_cand = sub_tbl.iloc[0]
        # Max POx
        max_POx.append(max_cand.P_Ox)
        # Best stuff
        best_offset.append(max_cand.separation)
        best_half.append(max_cand.half_light)
        best_mag.append(max_cand.ACS_i)
        best_ra.append(max_cand.ra)
        best_dec.append(max_cand.dec)
        best_Pc.append(max_cand.P_c)
        best_PO.append(max_cand.P_O)
        # PU
        PU.append(1.-np.sum(sub_tbl.P_Ox))
        iFRB.append(sub_tbl.iFRB.values[0])
    #
    stats = pandas.DataFrame(dict(max_POx=max_POx, best_offset=best_offset,
                                  best_half=best_half,
                                  best_mag=best_mag,
                                  best_ra=best_ra, best_dec=best_dec,
                                  best_PO=best_PO, best_Pc=best_Pc,
                                  PU=PU, iFRB=iFRB,
                                  ))
    stats.to_csv(stats_outfile)
    print("Wrote: {}".format(stats_outfile))

def analyze_box(truth_file, results_file, stats_file, figfile, secure_POx=0.9):
    # Load
    print("Loading..")
    #results = pandas.read_csv(results_file, index_col=0)
    stats = pandas.read_csv(stats_file, index_col=0)
    truth = pandas.read_csv(truth_file, index_col=0)

    # Merge
    merge = stats.join(truth)

    # Cuts
    secure = stats.max_POx > secure_POx

    # Max P(O|x)
    fig = plt.figure()
    with PdfPages(figfile) as pdf:


        # Histogram
        fig = plt.figure(figsize=(6,6))
        ax = plt.gca()
        _ = sns.histplot(stats, x='max_POx', ax=ax)
        ax.set_xlabel(r'$max[P(O|x)]$')
        ax.set_ylabel(r'PDF')

        pdf.savefig()
        plt.close()

        # Scatter P(O|x) truth
        fig = plt.figure(figsize=(6,6))
        _ = sns.displot(data=merge, x='mag_best', y='max_POx')#, ax=ax)
        pdf.savefig()
        plt.close()

        # Scatter P(O|x) assigned
        fig = plt.figure(figsize=(6,6))
        _ = sns.displot(data=merge, x='best_mag', y='max_POx')#, ax=ax)
        pdf.savefig()
        plt.close()

        # Scatter P_c, P(O|x) assigned
        fig = plt.figure(figsize=(6,6))
        _ = sns.displot(data=merge, x='best_Pc', y='max_POx')#, ax=ax)
        pdf.savefig()
        plt.close()

        # Scatter P_O, P(O|x) assigned
        fig = plt.figure(figsize=(6,6))
        _ = sns.displot(data=merge, x='best_PO', y='max_POx')#, ax=ax)
        pdf.savefig()
        plt.close()


        '''
        # Scatter mags
        fig = plt.figure(figsize=(6,6))
        ax = plt.gca()
        _ = sns.scatterplot(data=merge, x='mag_best', y='best_mag',
                            hue='max_POx', size='max_POx',
                            sizes=(1,50),
                            ax=ax)
        pdf.savefig()
        plt.close()
        '''

    print("Wrote: {}".format(figfile))


def run_it(prior, frb_ee,  frb_file, stats_file, outfile, stats_only=False,
           **kwargs):
    """
    Main run

    Returns:

    """
    # Run
    if not stats_only :
        run_bayes(frb_file,
              'SandBox/first_test_1000/cosmos_acs_iphot_200709.feather',
              prior, outfile, frb_ee=frb_ee, **kwargs)
    else:
        print("Running stats only")
    # Stats
    do_stats(outfile, stats_file)

def main(flg):

    # First 100k + adopted
    if flg & (2**0):
        # Prior
        prior = associate_defs.adopted.copy()
        prior['U'] = 0
        prior['O'] = 'inverse'
        # FRB
        frb_ee = {'a': 1., 'b': 1., 'theta': 0.}
        # Files
        sbox = 'mU_oU2_1arcsec'
        galaxies_file = os.path.join('SandBox', sbox, 'frbs_with_galaxies_100000.csv')
        frb_file = os.path.join('SandBox', sbox, 'frbs_100000.csv')
        results_file = os.path.join('SandBox', sbox, sbox+'_A_results.csv')
        stats_file = os.path.join('SandBox', sbox, sbox+'_A_stats.csv')

        # Run
        run_it(prior, frb_ee, frb_file, stats_file, results_file, stats_only=True)
        # Analyze
        analyze_box(galaxies_file, results_file, stats_file, sbox+'_A.pdf')

    # First 100k + conservative
    if flg & (2**1):
        # Prior
        prior = associate_defs.conservative.copy()
        prior['U'] = 0
        # FRB
        frb_ee = {'a': 1., 'b': 1., 'theta': 0.}
        # Files
        sbox = 'mU_oU2_1arcsec'
        galaxies_file = os.path.join('SandBox', sbox, 'frbs_with_galaxies_100000.csv')
        frb_file = os.path.join('SandBox', sbox, 'frbs_100000.csv')
        results_file = os.path.join('SandBox', sbox, sbox+'_C_results.csv')
        stats_file = os.path.join('SandBox', sbox, sbox+'_C_stats.csv')

        # Run
        run_it(prior, frb_ee, frb_file, stats_file, results_file, stats_only=False)
        # Analyze
        analyze_box(galaxies_file, results_file, stats_file, sbox+'_C.pdf')

    # First 100k + adopted + U0.05
    if flg & (2**2):
        # Prior
        prior = associate_defs.adopted.copy()
        prior['U'] = 0.05
        prior['O'] = 'inverse'
        # FRB
        frb_ee = {'a': 1., 'b': 1., 'theta': 0.}
        # Files
        sbox = 'mU_oU2_1arcsec'
        galaxies_file = os.path.join('SandBox', sbox, 'frbs_with_galaxies_100000.csv')
        frb_file = os.path.join('SandBox', sbox, 'frbs_100000.csv')
        results_file = os.path.join('SandBox', sbox, sbox+'_A0.05_results.csv')
        stats_file = os.path.join('SandBox', sbox, sbox+'_A0.05_stats.csv')

        # Run
        run_it(prior, frb_ee, frb_file, stats_file, results_file, stats_only=False)
        # Analyze
        analyze_box(galaxies_file, results_file, stats_file, sbox+'_A0.05.pdf')

    # SB-2, 3, 4
    if flg & (2**3):
        # Prior
        iprior = associate_defs.adopted.copy()

        # Make em
        prior1 = iprior.copy()
        prior1['theta']['method'] = 'exp'
        prior2 = iprior.copy()
        prior2['theta']['method'] = 'core'
        prior3 = iprior.copy()
        prior3['theta']['method'] = 'uniform'
        prior3['theta']['max'] = 2.

        # Files
        for ss, prior, sbox in zip(
                np.arange(3),
                [prior1, prior2, prior3],
                ['mag_20-23_oExp_locU_0.1-1', 'mag_20-23_oCore_locU_0.1-1',
                     'mag_20-23_oU2_locU_0.1-1']):
            #
            if ss < 2:
                continue
            #
            galaxies_file = glob.glob(os.path.join('SandBox', sbox, 'frbs_with_galaxies*'))[0]
            frb_file = glob.glob(os.path.join('SandBox', sbox, 'frbs_46*'))[0]
            #
            results_file = os.path.join('SandBox', sbox, sbox+'_A_results.csv')
            stats_file = os.path.join('SandBox', sbox, sbox+'_A_stats.csv')

            # Run
            run_it(prior, 'loc_sig', frb_file, stats_file, results_file, stats_only=False)
            # Analyze
            analyze_box(galaxies_file, results_file, stats_file, sbox+'_A.pdf')


    # PU10
    if flg & (2**4):
        # Prior
        prior = associate_defs.adopted.copy()
        prior['U'] = 0.10
        # FRB
        # Files
        sbox = 'PU10'

        # m<23
        galaxies_file = os.path.join('SandBox', sbox, 'frbs_with_galaxies_50000_mag_20-23_oU2_locU_0.1-1_PU10.csv')
        frb_file = os.path.join('SandBox', sbox, 'frbs_50000_mag_20-23_oU2_locU_0.1-1_PU10.csv')
        results_file = os.path.join('SandBox', sbox, sbox+'mag23_A0.10_results.csv')
        stats_file = os.path.join('SandBox', sbox, sbox+'mag23_A0.10_stats.csv')

        # Run
        run_it(prior, 'loc_sig', frb_file, stats_file, results_file, stats_only=True)
        # Analyze
        analyze_box(galaxies_file, results_file, stats_file, sbox+'mag23_A0.10.pdf')

        # m<25
        galaxies_file = os.path.join('SandBox', sbox, 'frbs_with_galaxies_50000_mag_20-25_oU2_locU_0.1-1_PU10.csv')
        frb_file = os.path.join('SandBox', sbox, 'frbs_50000_mag_20-25_oU2_locU_0.1-1_PU10.csv')
        results_file = os.path.join('SandBox', sbox, sbox + 'mag25_A0.10_results.csv')
        stats_file = os.path.join('SandBox', sbox, sbox + 'mag25_A0.10_stats.csv')

        # Run
        run_it(prior, 'loc_sig', frb_file, stats_file, results_file, stats_only=True)
        # Analyze
        analyze_box(galaxies_file, results_file, stats_file, sbox + 'mag25_A0.10.pdf')

    # First 100k + adopted with inverse1
    if flg & (2**5):
        # Prior
        prior = associate_defs.adopted.copy()
        prior['U'] = 0
        prior['O'] = 'inverse1'
        # FRB
        frb_ee = {'a': 1., 'b': 1., 'theta': 0.}
        # Files
        sbox = 'mU_oU2_1arcsec'
        galaxies_file = os.path.join('SandBox', sbox, 'frbs_with_galaxies_100000.csv')
        frb_file = os.path.join('SandBox', sbox, 'frbs_100000.csv')
        # Outputs
        results_file = os.path.join('SandBox', sbox, sbox+'_A1_results.csv')
        stats_file = os.path.join('SandBox', sbox, sbox+'_A1_stats.csv')

        # Run
        run_it(prior, frb_ee, frb_file, stats_file, results_file, stats_only=False)
        # Analyze
        analyze_box(galaxies_file, results_file, stats_file, sbox+'_A1.pdf')

    # First 100k + adopted with inverse2
    if flg & (2**6):
        # Prior
        prior = associate_defs.adopted.copy()
        prior['U'] = 0
        prior['O'] = 'inverse2'
        # FRB
        frb_ee = {'a': 1., 'b': 1., 'theta': 0.}
        # Files
        sbox = 'mU_oU2_1arcsec'
        galaxies_file = os.path.join('SandBox', sbox, 'frbs_with_galaxies_100000.csv')
        frb_file = os.path.join('SandBox', sbox, 'frbs_100000.csv')
        # Outputs
        results_file = os.path.join('SandBox', sbox, sbox+'_A2_results.csv')
        stats_file = os.path.join('SandBox', sbox, sbox+'_A2_stats.csv')

        # Run
        run_it(prior, frb_ee, frb_file, stats_file, results_file, stats_only=False)
        # Analyze
        analyze_box(galaxies_file, results_file, stats_file, sbox+'_A2.pdf')

    # SB-1 + inverse + uniform offset w/ theta_max = 2phi
    if flg & (2**7):
        # Prior
        prior = associate_defs.conservative.copy()
        prior['U'] = 0
        prior['O'] = 'inverse'
        prior['theta']['max'] = 2. # phi
        # FRB
        frb_ee = {'a': 1., 'b': 1., 'theta': 0.}
        # Files
        sbox = 'mU_oU2_1arcsec'
        galaxies_file = os.path.join('SandBox', sbox, 'frbs_with_galaxies_100000.csv')
        frb_file = os.path.join('SandBox', sbox, 'frbs_100000.csv')
        # Outputs
        results_file = os.path.join('SandBox', sbox, sbox+'_M_results.csv')
        stats_file = os.path.join('SandBox', sbox, sbox+'_M_stats.csv')

        # Run
        run_it(prior, frb_ee, frb_file, stats_file, results_file, stats_only=False)
        # Analyze
        analyze_box(galaxies_file, results_file, stats_file, sbox+'_M.pdf')

    # Unseen
    if flg & (2**8):
        # Prior
        prior = associate_defs.adopted.copy()
        prior['U'] = 0
        prior['O'] = 'inverse'
        # FRB
        frb_ee = {'a': 1., 'b': 1., 'theta': 0.}
        # Files
        sbox = 'unseen'
        galaxies_file = os.path.join('SandBox', sbox, 'frbs_with_galaxies_100000.csv')
        frb_file = os.path.join('SandBox', sbox, 'frbs_100000.csv')

        # Outputs
        results_file = os.path.join('SandBox', sbox, sbox+'_AU0_results.csv')
        stats_file = os.path.join('SandBox', sbox, sbox+'_AU0_stats.csv')
        # Run
        run_it(prior, frb_ee, frb_file, stats_file, results_file, stats_only=False,
               mag_lim=25.)
        # Analyze
        analyze_box(galaxies_file, results_file, stats_file, sbox+'_AU0.pdf')

        prior['U'] = 0.1
        # Outputs
        results_file = os.path.join('SandBox', sbox, sbox+'_AU0.1_results.csv')
        stats_file = os.path.join('SandBox', sbox, sbox+'_AU0.1_stats.csv')
        # Run
        run_it(prior, frb_ee, frb_file, stats_file, results_file, stats_only=False,
               mag_lim=25.)
        # Analyze
        analyze_box(galaxies_file, results_file, stats_file, sbox+'_AU0.1.pdf')



# Command line execution
if __name__ == '__main__':

    if len(sys.argv) == 1:
        flg = 0
        flg += 2**0   # First 100k + inverse
        #flg += 2**1   # First 100k + conservative
        #flg += 2**2   # First 100k + adopted + P(U)=0.05
        #flg += 2**3   # SB-2, 3, 4
        #flg += 2**4   # PU10
        #flg += 2**5   # First 100k + inverse1
        #flg += 2**6   # First 100k + inverse2
        #flg += 2**7   # SB-1 match with inverse
        #flg += 2**8   # Unseen FRBs

    else:
        flg = int(sys.argv[1])

    #
    main(flg)
