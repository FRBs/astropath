
import numpy as np
import sys

import matplotlib as mpl
import seaborn as sns

import pandas

#mpl.rcParams['font.family'] = 'stixgeneral'

from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec


from statsmodels.stats.proportion import proportion_confint


from astropy import units
from astropy.coordinates import SkyCoord

from IPython import embed

# Local
#sys.path.append(os.path.abspath("../Analysis/py"))
#import associate_defs
#import analysis

# Globals
handletextpad=0.3

# Variables
cpoffset = r'$P(\omega|O)$'
cPOi = r'$P(O_i)$'
cPOxi = r'$P(O_i|x)$'
chalf = r'$\phi_{1/2}$'
cmhalf = '\phi_{1/2}'


# COSMOS
plate_scale = 0.05

def fig_sb_posteriors(outfile='fig_sb_posteriors.png'):
    """

    Args:
        outfile:

    Returns:

    """
    set_mplrc()

    sns.set_theme()
    sns.set_style('whitegrid')
    sns.set_context('paper')

    # Load
    print("Loading...")
    truth_file = '../Analysis/SandBox/first_test_1000/frbs_100k_with_galaxies.csv'
    stats_file = '../Analysis/SandBox/first_test_1000/first_test_100k_U_stats.csv'
    results_file = '../Analysis/SandBox/first_test_1000/first_test_100k_U_results.csv'

    stats = pandas.read_csv(stats_file, index_col=0)
    truth = pandas.read_csv(truth_file, index_col=0)
    results = pandas.read_csv(results_file, index_col=0)


    # Merge
    merge = stats.join(truth)

    # Plot
    plt.figure(figsize=(6, 12))
    gs = gridspec.GridSpec(2,1)

    ax = plt.subplot(gs[0])
    above_min = results[results.P_Ox > 0.01]
    print("Number of candidates with 0 > P_Ox < 0.01 = {} ".format(
        np.sum( (results.P_Ox < 0.01) & (results.P_Ox > 0))))
    print("Of {} total candidates".format(len(results)))
    _ = sns.histplot(above_min, x='P_Ox', ax=ax, stat='probability', color='green')
    print(np.sum((results.P_Ox < 0.01) & (results.P_Ox > 0))/len(results))
    ax.set_xlabel(r'$P(O_i|x)$')
    ax.set_ylabel(r'PDF')
    ax.minorticks_on()

    ylbl = 0.9
    lsz = 21.
    ax.text(0.50, ylbl, '(a)', transform=ax.transAxes, size=lsz, ha='center')

    # Max P(O|x)
    ax1 = plt.subplot(gs[1])
    _ = sns.histplot(stats, x='max_POx', ax=ax1, stat='probability')
    ax1.set_xlabel(r'$max[P(O_i|x)]$')
    ax1.set_ylabel(r'PDF')
    ax1.minorticks_on()

    ax1.text(0.50, ylbl, '(b)', transform=ax1.transAxes, size=lsz, ha='center')

    # Stats
    psecure = 0.95
    print("Fraction of FRBs with max[P(O|x)] > {}: {}".format(psecure,
        np.sum(stats.max_POx > psecure)/len(stats)))

    print("Fraction of FRBs with P(O|x) > 0.90 or <0.10: {}".format(
        np.sum(np.any([results.P_Ox.values > 0.90, results.P_Ox.values < 0.1],axis=0))/len(results)))

    # Values
    '''
    ax.text(0.55, 0.80, r'$r_0 = {:0.2f}'.format(r0_ML2)+'_{-'+'{:0.2f}'.format(
        r0_ML2-rerr[0])+'}^{+'+'{:0.2f}'.format(rerr[1]-r0_ML2)+'} \; (h_{100}^{-1} \, $ Mpc)',
            transform=ax.transAxes, size='large', ha='left')
    '''

    # Label me
    #ax.xaxis.set_major_locator(plt.MultipleLocator(1.))
    #ax.xaxis.set_major_formatter(FormatStrFormatter(r'$%.3f$'))
    #ax.set_ylim(0., 1)
    #ax.set_xlim(0., 1.)


    for iax in [ax,ax1]:
        set_fontsize(iax, 15.)

    # End
    plt.tight_layout(pad=0.3, h_pad=0.2, w_pad=0.2)
    print('Writing {:s}'.format(outfile))
    kwargs = {}
    if 'png' in outfile:
        kwargs['dpi'] = 500
    plt.savefig(outfile, **kwargs)
    plt.close()


def fig_sb_PP(outfile='fig_sb_PP.png', exploring=False, SB3=False):
    """

    Args:
        outfile:

    Returns:

    """
    if exploring:
        outfile='fig_sb_PP_exploring.png'
    elif SB3:
        outfile = 'fig_sb_PP_ASKAP.png'
    set_mplrc()

    # Load
    print("Loading...")

    if SB3:
        truth_file = '../Analysis/SandBox/mag_20-23_oExp_locU_0.1-1/frbs_with_galaxies_46699_mag_20-23_oExp_locU_0.1-1.csv'
    else:
        truth_file = '../Analysis/SandBox/mU_oU2_1arcsec/frbs_with_galaxies_100000.csv'
    truth = pandas.read_csv(truth_file, index_col=0)

    if exploring:
        stats_fileU = '../Analysis/SandBox/mU_oU2_1arcsec/mU_oU2_1arcsec_A_stats.csv'
        statsU = pandas.read_csv(stats_fileU, index_col=0)

        stats_fileU05 = '../Analysis/SandBox/mU_oU2_1arcsec/mU_oU2_1arcsec_A1_stats.csv'
        statsU05 = pandas.read_csv(stats_fileU05, index_col=0)

        stats_fileC = '../Analysis/SandBox/mU_oU2_1arcsec/mU_oU2_1arcsec_A2_stats.csv'
        statsC = pandas.read_csv(stats_fileC, index_col=0)

        #stats_fileM = '../Analysis/SandBox/mU_oU2_1arcsec/mU_oU2_1arcsec_M_stats.csv'
        #statsM = pandas.read_csv(stats_fileM, index_col=0)
    elif SB3:
        stats_fileU = '../Analysis/SandBox/mag_20-23_oExp_locU_0.1-1/mag_20-23_oExp_locU_0.1-1_A2_stats.csv'
        statsU = pandas.read_csv(stats_fileU, index_col=0)
    else:
        stats_fileU = '../Analysis/SandBox/mU_oU2_1arcsec/mU_oU2_1arcsec_A2_stats.csv'
        statsU = pandas.read_csv(stats_fileU, index_col=0)

        stats_fileU05 = '../Analysis/SandBox/mU_oU2_1arcsec/mU_oU2_1arcsec_A20.05_stats.csv'
        statsU05 = pandas.read_csv(stats_fileU05, index_col=0)

        stats_fileC = '../Analysis/SandBox/mU_oU2_1arcsec/mU_oU2_1arcsec_C_stats.csv'
        statsC = pandas.read_csv(stats_fileC, index_col=0)

    # Coords
    true_coords = SkyCoord(ra=truth.ra.values, dec=truth.dec.values, unit='deg')

    pred_coordsU = SkyCoord(ra=statsU.best_ra.values, dec=statsU.best_dec.values, unit='deg')
    if not SB3:
        pred_coordsU05 = SkyCoord(ra=statsU05.best_ra.values, dec=statsU05.best_dec.values, unit='deg')
        pred_coordsC = SkyCoord(ra=statsC.best_ra.values, dec=statsC.best_dec.values, unit='deg')

    # Match?
    if SB3:
        pred_list = [pred_coordsU]
        stat_list = [statsU]
    elif exploring:
        #pred_coordsM = SkyCoord(ra=statsM.best_ra.values, dec=statsM.best_dec.values, unit='deg')
        pred_list = [pred_coordsU, pred_coordsU05, pred_coordsC]#, pred_coordsM]
        stat_list = [statsU, statsU05, statsC]#, statsM]
    else:
        pred_list = [pred_coordsU, pred_coordsU05, pred_coordsC]
        stat_list = [statsU, statsU05, statsC]
    for pred, stats in zip(pred_list, stat_list):
        sep = pred.separation(true_coords).to('arcsec')
        correct = sep < 0.01*units.arcsec
        stats['correct'] = correct

    print(np.max(stats.max_POx))
    nbin = 10
    nsub = len(truth) // nbin
    if exploring:
        dictU, dictU05, dictC, dictM = dict(lbl=r'$\Sigma(m)$'), \
                                dict(lbl=r'$\Sigma(m) \phi$'), \
                                       dict(lbl=r'$\Sigma(m) \phi^2$'), \
                                            dict(lbl=r'$U + \Sigma(m)$')
    elif SB3:
        dictU = dict(lbl='Adopted')
    else:
        dictU, dictU05, dictC = dict(lbl='Adopted'), \
                                dict(lbl=r'Adopted, $P(U)=0.05$'), \
                                dict(lbl='Conservative')

    # Stats time
    if SB3:
        dict_list = [dictU]
    elif exploring:
        dict_list = [dictU, dictU05, dictC]#, dictM]
    else:
        dict_list = [dictU, dictU05, dictC]
    for pdict, stats in zip(dict_list, stat_list):
        # Init dict
        pdict['x0'] = np.zeros(nbin)
        pdict['x1'] = np.zeros(nbin)
        pdict['perc'] = np.zeros(nbin)
        pdict['ci_low'] = np.zeros(nbin)
        pdict['ci_upp'] = np.zeros(nbin)
        # Sort
        isrt = np.argsort(stats.max_POx.values)
        for ss in range(nbin):
            x0 = stats.max_POx[isrt[ss*nsub]]
            x1 = stats.max_POx[isrt[(ss+1)*nsub-1]]
            pdict['x0'][ss] = x0
            pdict['x1'][ss] = x1
            #
            in_bin = (stats.max_POx > x0) & (stats.max_POx <= x1)
            nobj = int(np.sum(in_bin))
            nhits = int(np.sum(stats.correct.values[in_bin]))
            try:
                pdict['perc'][ss] = nhits/nobj
            except:
                embed(header='211 of figs')
            # Stats
            pdict['ci_low'][ss], pdict['ci_upp'][ss] = proportion_confint(nhits, nobj)

    # Plot
    plt.figure(figsize=(6, 5))
    gs = gridspec.GridSpec(1,1)
    ax = plt.subplot(gs[0])

    # 1-1
    ax.plot([0.,1.], [0., 1.], 'k--', label='1-1')

    # Models
    for pdict in dict_list:
        xval = (pdict['x0']+pdict['x1'])/2.
        ax.errorbar(xval, pdict['perc'],
                    xerr=xval-pdict['x0'],
                    #yerr=[pdict['perc']-pdict['ci_low'],
                    #      pdict['ci_upp']-pdict['perc']],
                    fmt='.',
                    label=pdict['lbl'], capthick=2, capsize=2)

    # Label me
    ax.set_xlabel(r'max $P(O_i|X)$')
    ax.set_ylabel('Fraction correct')
    ax.minorticks_on()
    #ax.xaxis.set_major_locator(plt.MultipleLocator(1.))
    #ax.xaxis.set_major_formatter(FormatStrFormatter(r'$%.3f$'))
    ax.set_ylim(0., 1.05)
    ax.set_xlim(0., 1.05)

    # Legend
    legend = ax.legend(loc='lower right', scatterpoints=1, borderpad=0.2,
                       handletextpad=handletextpad, fontsize='large')

    set_fontsize(ax, 15.)
    # End
    plt.tight_layout(pad=0.2, h_pad=0., w_pad=0.1)
    print('Writing {:s}'.format(outfile))
    kwargs = {}
    if 'png' in outfile:
        kwargs['dpi'] = 500
    plt.savefig(outfile, **kwargs)
    plt.close()


def fig_sb_PU(outfile='fig_sb_PU.png', exploring=False, SB3=False, mag=23):
    """

    Args:
        outfile:

    Returns:

    """
    set_mplrc()

    # Load
    print("Loading...")

    #truth23_file = '../Analysis/SandBox/PU10/frbs_with_galaxies_50000_mag_20-23_oU2_locU_0.1-1_PU10.csv'
    stats23_file = '../Analysis/SandBox/PU10/PU10mag23_A0.10_stats.csv'
    #truth25_file = '../Analysis/SandBox/PU10/frbs_with_galaxies_50000_mag_20-25_oU2_locU_0.1-1_PU10.csv'
    stats25_file = '../Analysis/SandBox/PU10/PU10mag25_A0.10_stats.csv'

    stats23 = pandas.read_csv(stats23_file, index_col=0)
    stats25 = pandas.read_csv(stats25_file, index_col=0)
    stat_list = [stats23, stats25]

    # Coords
    dict23, dict25  = dict(lbl='mag23'), \
                            dict(lbl='mag25')
    dict_list = [dict23, dict25]

    binv = [0., 0.05, 0.1, 0.3, 0.9, 1.]
    nbin = len(binv)-1
    nsub = len(stats23) // nbin
    for pdict, stats in zip(dict_list, stat_list):
        # Init dict
        pdict['x0'] = np.zeros(nbin)
        pdict['x1'] = np.zeros(nbin)
        pdict['perc'] = np.zeros(nbin)
        pdict['ci_low'] = np.zeros(nbin)
        pdict['ci_upp'] = np.zeros(nbin)
        # Sort
        isrt = np.argsort(stats.PU.values)
        for ss in range(nbin):
            x0 = binv[ss]
            x1 = binv[ss+1]
            pdict['x0'][ss] = x0
            pdict['x1'][ss] = x1
            #
            in_bin = (stats.PU > x0) & (stats.PU <= x1)
            print("x0,x1,Nbin: {},{},{}".format(x0, x1, np.sum(in_bin)))
            # Last 5000 are unseen
            nobj = int(np.sum(in_bin))
            nU = int(np.sum(stats[in_bin].iFRB > 45000))
            pdict['perc'][ss] = nU/nobj
            # Stats
            pdict['ci_low'][ss], pdict['ci_upp'][ss] = proportion_confint(nU, nobj)

    # Plot
    plt.figure(figsize=(6, 5))
    gs = gridspec.GridSpec(1,1)
    ax = plt.subplot(gs[0])

    # 1-1
    ax.plot([0.,1.], [0., 1.], 'k--', label='1-1')

    # Models
    for pdict in dict_list:
        xval = (pdict['x0']+pdict['x1'])/2.
        ax.errorbar(xval, pdict['perc'],
                xerr=xval-pdict['x0'],
                    #yerr=[pdict['perc']-pdict['ci_low'],
                    #      pdict['ci_upp']-pdict['perc']],
                fmt='.',
                label=pdict['lbl'], capthick=2, capsize=2)

    # Label me
    ax.set_xlabel(r'$P(U)$')
    ax.set_ylabel('Fraction Unseen')
    ax.minorticks_on()
    #ax.xaxis.set_major_locator(plt.MultipleLocator(1.))
    #ax.xaxis.set_major_formatter(FormatStrFormatter(r'$%.3f$'))
    ax.set_ylim(0., 1.05)
    ax.set_xlim(0., 1.05)

    # Legend
    legend = ax.legend(loc='lower right', scatterpoints=1, borderpad=0.2,
                       handletextpad=handletextpad, fontsize='large')

    set_fontsize(ax, 15.)
    # End
    plt.tight_layout(pad=0.2, h_pad=0., w_pad=0.1)
    print('Writing {:s}'.format(outfile))
    kwargs = {}
    if 'png' in outfile:
        kwargs['dpi'] = 500
    plt.savefig(outfile, **kwargs)
    plt.close()



def fig_sb_COSMOS(outfile='fig_sb_COSMOS.png'):
    """

    Args:
        outfile:

    Returns:

    """
    sns.set_theme()
    sns.set_style('whitegrid')
    sns.set_context('paper')

    # Load
    print("Loading...")
    galaxy_cat_file = '../Analysis/SandBox/first_test_1000/cosmos_acs_iphot_200709.feather'
    sb_galaxies = pandas.read_feather(galaxy_cat_file)
    # Cut
    sb_galaxies = sb_galaxies[np.isfinite(sb_galaxies.a_image)]

    # Scale
    sb_galaxies.a_image *= plate_scale


    # Plot
    plt.figure(figsize=(8, 8))
    ax = plt.gca()

    _ = sns.histplot(data=sb_galaxies, x='mag_best', y='a_image', ax=ax,
                     log_scale=(False,True))

    ax.set_xlabel(r'$m_r$  (AB mag)')
    ax.set_ylabel(r'$\phi$ (arcsec)')
    ax.minorticks_on()

    # Values
    '''
    ax.text(0.55, 0.80, r'$r_0 = {:0.2f}'.format(r0_ML2)+'_{-'+'{:0.2f}'.format(
        r0_ML2-rerr[0])+'}^{+'+'{:0.2f}'.format(rerr[1]-r0_ML2)+'} \; (h_{100}^{-1} \, $ Mpc)',
            transform=ax.transAxes, size='large', ha='left')
    '''

    # Label me
    ax.xaxis.set_major_locator(plt.MultipleLocator(5.))
    #ax.set_ylim(0., 15.)
    ax.set_xlim(10., 32.)

    set_fontsize(ax, 17.)

    # End
    #plt.tight_layout(pad=0.3, h_pad=0., w_pad=0.2)
    print('Writing {:s}'.format(outfile))
    kwargs = {}
    if 'png' in outfile:
        kwargs['dpi'] = 500
    plt.savefig(outfile, **kwargs)
    plt.close()

def fig_sb_mag(outfile='fig_sb_POx_mags.png', Psecure=0.95):

    # Load up
    truth_file = '../Analysis/SandBox/mU_oU2_1arcsec/frbs_with_galaxies_100000.csv'
    truth = pandas.read_csv(truth_file, index_col=0)

    true_coords = SkyCoord(ra=truth.ra.values, dec=truth.dec.values, unit='deg')

    stats_file = '../Analysis/SandBox/mU_oU2_1arcsec/mU_oU2_1arcsec_A_stats.csv'
    stats = pandas.read_csv(stats_file, index_col=0)

    # Correct
    pred_coords = SkyCoord(ra=stats.best_ra.values, dec=stats.best_dec.values, unit='deg')
    sep = pred_coords.separation(true_coords).to('arcsec')
    correct = sep < 0.01 * units.arcsec
    stats['correct'] = correct

    # Mag cuts
    mag_cuts = [21., 23., 24, 25., 26.]
    ncuts = len(mag_cuts)-1

    # Plot
    plt.figure(figsize=(6, 5))
    gs = gridspec.GridSpec(2,2)

    for cut in range(ncuts):
        ax = plt.subplot(gs[cut])

        # Cut
        mag_cut = (truth.mag_best <= mag_cuts[cut+1]) & (truth.mag_best > mag_cuts[cut])
        # Stats
        secure = stats[mag_cut].max_POx > Psecure
        TP = secure & stats[mag_cut].correct

        ax = sns.histplot(data=stats[mag_cut], x='max_POx', hue='correct')
        #
        lsz = 9.
        ax.text(0.1, 0.3, 'TP({:0.2f}) = {:0.2f}'.format(Psecure, np.sum(TP) / np.sum(secure)),
                transform=ax.transAxes, size=lsz, ha='left')
        ax.text(0.1, 0.4, 'f_secure({:0.2f}) = {:0.3f}'.format(Psecure, np.sum(secure) / np.sum(mag_cut)),
                transform=ax.transAxes, size=lsz, ha='left')
        #
        ax.set_yscale('log')

        set_fontsize(ax, 11.)
        ax.set_title('mag=[{},{}]'.format(mag_cuts[cut], mag_cuts[cut+1]))

    # End
    plt.tight_layout(pad=0.2, h_pad=0., w_pad=0.1)
    print('Writing {:s}'.format(outfile))
    kwargs = {}
    if 'png' in outfile:
        kwargs['dpi'] = 500
    plt.savefig(outfile, **kwargs)
    plt.close()


def read_unseen(name, truth='frbs_with_galaxies_50000_mag_20-25_oU2_locU_0.1-1_PU10.csv'):
    # read stats and results of SB with unseen galaxies and return results only for unseen
    # FRB hosts
    truth_file = truth
    stats_file = f'{name}_stats.csv'
    results_file = f'{name}_results.csv'

    stats = pandas.read_csv(stats_file, index_col=0)
    truth = pandas.read_csv(truth_file, index_col=0)
    results = pandas.read_csv(results_file, index_col=0, keep_default_na=False)
    unseen_results = results[results['iFRB'] > truth.index.max()]
    unseen_stats = stats[stats['iFRB'] > truth.index.max()]
    return unseen_stats, unseen_results, truth


def fig_unseen_hosts_cdf(outfile='fig_unseen_hosts_cdf.png'):
    unseen_stats_C, unseen_results_C, truth_file = read_unseen(name='PU10mag25_C')
    unseen_stats_A, unseen_results_A, truth_file = read_unseen(name='PU10mag25_A0')

    psecure = 0.95

    print(f"Fraction of secure for Conservative (for unseen host FRBs): "
          f"{(unseen_stats_C['max_POx']  > psecure).sum()/len(unseen_stats_C)}")

    print(f"Fraction of secure for Adopted (for unseen host FRBs): "
          f"{(unseen_stats_A['max_POx']  > psecure).sum()/len(unseen_stats_A)}")

    unseen_stats_C['theta_over_phi'] = unseen_stats_C['best_offset'] / unseen_stats_C['best_half']
    unseen_stats_A['theta_over_phi'] = unseen_stats_A['best_offset'] / unseen_stats_A['best_half']

    secure_prob = 0.95
    unseen_stats_C_secure = unseen_stats_C[unseen_stats_C['max_POx'] > secure_prob]
    unseen_stats_A_secure = unseen_stats_A[unseen_stats_A['max_POx'] > secure_prob]

    hist_C, bins = np.histogram(unseen_stats_C['theta_over_phi'], bins=500)
    hist_C_secure, bins = np.histogram(unseen_stats_C_secure['theta_over_phi'], bins=bins)

    cdf_C = np.cumsum(hist_C)/len(unseen_stats_C)
    cdf_C_secure = np.cumsum(hist_C_secure)/len(unseen_stats_C_secure)

    hist_A, bins = np.histogram(unseen_stats_A['theta_over_phi'], bins=bins)
    hist_A_secure, bins = np.histogram(unseen_stats_A_secure['theta_over_phi'], bins=bins)

    cdf_A = np.cumsum(hist_A)/len(unseen_stats_A)
    cdf_A_secure = np.cumsum(hist_A_secure)/len(unseen_stats_A_secure)

    fig, ax1 = plt.subplots(1,1)
    ax1.plot(bins[:-1], cdf_C, color='g', label='Conservative')
    ax1.plot(bins[:-1], cdf_C_secure, color='g', linestyle='--')

    ax1.plot(bins[:-1], cdf_A, color='b', label='Adopted')
    ax1.plot(bins[:-1], cdf_A_secure, color='b', linestyle='--')

    ax1.set_xlim([0, 6])
    ax1.set_ylim([0, 0.13])
    ax1.set_xlabel(r'$\theta/\phi$')
    ax1.set_ylabel(r'CDF')
    ax1.minorticks_on()
    plt.grid()
    plt.legend()
    plt.savefig(outfile)


def set_fontsize(ax,fsz):
    '''
    Parameters
    ----------
    ax : Matplotlib ax class
    fsz : float
      Font size
    '''
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(fsz)

def set_mplrc():
    mpl.rcParams['mathtext.default'] = 'it'
    mpl.rcParams['font.size'] = 12
    #mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    #mpl.rc('font',family='Times New Roman')
    #mpl.rcParams['text.latex.preamble'] = [r'\boldmath']
    mpl.rcParams['mathtext.fontset'] = 'custom'
    mpl.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
    mpl.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
    mpl.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
    mpl.rc('text', usetex=True)


#### ########################## #########################
def main(flg_fig):

    if flg_fig == 'all':
        flg_fig = np.sum( np.array( [2**ii for ii in range(25)] ))
    else:
        flg_fig = int(flg_fig)

    # P(theta)
    if flg_fig & (2**0):
        fig_sb_posteriors()

    # PP - SB-1
    if flg_fig & (2**1):
        fig_sb_PP()

    # COSMOS
    if flg_fig & (2**2):
        fig_sb_COSMOS()

    # PP for the prior
    if flg_fig & (2**3):
        fig_sb_PP(exploring=True)

    # PP - SB-3
    if flg_fig & (2**4):
        fig_sb_PP(SB3=True)

    # PP - SB-U
    if flg_fig & (2**5):
        fig_sb_PU()

    # P(O|x) for SB-1 and mag cuts
    if flg_fig & (2**6):
        fig_sb_mag()


# Command line execution
if __name__ == '__main__':


    if len(sys.argv) == 1:
        flg_fig = 0
        #flg_fig += 2**0   # Posteriors
        #flg_fig += 2**1   # PP
        #flg_fig += 2**2   # COSMOS
        #flg_fig += 2**3   # Figuring out the 'best' prior
        #flg_fig += 2**4   # PP ASKAP
        #flg_fig += 2**5   # PU
        flg_fig += 2**6   # Stats with magnitude
    else:
        flg_fig = sys.argv[1]

    main(flg_fig)