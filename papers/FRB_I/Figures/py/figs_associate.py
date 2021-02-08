
import numpy as np
import glob, os, sys, json
import pdb

#import healpy as hp

from cycler import cycler
import matplotlib as mpl
import seaborn as sns

import pandas

from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

from scipy.interpolate import interp1d
from scipy.stats import kstest
from scipy.signal import convolve


from astropy import units

from photutils import SkyEllipticalAperture

from frb.associate import frbassociate
from astropath import bayesian
from frb.associate import frbs
from frb.dm import igm

from IPython import embed

# Local
sys.path.append(os.path.abspath("../Analysis/py"))
import associate_defs
import analysis

# Globals
handletextpad=0.3

# Variables
cpoffset = r'$p(\omega|O)$'
cPOi = r'$P(O_i)$'
cPOxi = r'$P(O_i|x)$'
chalf = r'$\phi_{1/2}$'
#cmhalf = '\phi_{1/2}'
cmhalf = '\phi'



def fig_theta_priors(outfile='fig_theta_priors.png'):
    """

    Args:
        outfile:

    Returns:

    """
    set_mplrc()

    theta_max = associate_defs.theta_max
    theta_u = dict(method='uniform', max=theta_max)
    theta_c = dict(method='core', max=theta_max)
    theta_e = dict(method='exp', max=theta_max)

    thetas = np.linspace(0., theta_max, 1000)
    dtheta = thetas[1] - thetas[0]

    # Plot
    plt.figure(figsize=(6, 5))
    gs = gridspec.GridSpec(1,1)
    ax = plt.subplot(gs[0])

    # Uniform
    ptheta_u = bayesian.pw_Oi(thetas, 1., theta_u)
    ax.plot(thetas, ptheta_u, label='Uniform', color='k')
    norm = theta_max

    # rcore
    ptheta_c = bayesian.pw_Oi(thetas, 1., theta_c)
    scl_c = norm / np.sum(ptheta_c*dtheta)
    scl_c = 1.
    #
    ax.plot(thetas, scl_c*ptheta_c, label='Core', color='b')

    # exponential
    ptheta_e = bayesian.pw_Oi(thetas, 1., theta_e)
    scl_e = norm / np.sum(ptheta_e*dtheta)
    scl_e = 1.
    #
    ax.plot(thetas, scl_e*ptheta_e, label='Exponential', color='g')


    # Values
    '''
    ax.text(0.55, 0.80, r'$r_0 = {:0.2f}'.format(r0_ML2)+'_{-'+'{:0.2f}'.format(
        r0_ML2-rerr[0])+'}^{+'+'{:0.2f}'.format(rerr[1]-r0_ML2)+'} \; (h_{100}^{-1} \, $ Mpc)',
            transform=ax.transAxes, size='large', ha='left')
    '''

    # Label me
    ax.set_xlabel(r'$\theta/'+cmhalf+'$')
    ax.set_ylabel(cpoffset)
    ax.xaxis.set_major_locator(plt.MultipleLocator(1.))
    ax.minorticks_on()
    #ax.xaxis.set_major_formatter(FormatStrFormatter(r'$%.3f$'))
    #ax.set_ylim(0., 0.1)

    # Legend
    legend = ax.legend(loc='upper right', scatterpoints=1, borderpad=0.2,
                       handletextpad=handletextpad, fontsize='large')

    set_fontsize(ax, 15.)
    # End
    plt.tight_layout(pad=0.2, h_pad=0., w_pad=0.1)
    print('Writing {:s}'.format(outfile))
    kwargs = {}
    if 'png' in outfile:
        kwargs['dpi'] = 700
    plt.savefig(outfile, **kwargs)
    plt.close()


def fig_prior_vs_posterior(outfile='fig_prior_vs_posterior.png',
                           POx_secure=associate_defs.POx_secure,
                           prior_mode='conservative'):
    """

    Args:
        outfile:

    Returns:

    """
    set_mplrc()

    frb_pre = analysis.get_candidates()

    # Init prior
    prior = getattr(associate_defs, prior_mode).copy()

    thetas = np.linspace(0., associate_defs.theta_max, 1000)
    dtheta = thetas[1] - thetas[0]

    bins_theta = np.linspace(0., 6., 20)

    # Plot
    plt.figure(figsize=(6, 5))
    gs = gridspec.GridSpec(2,2)


    # Loop on theta priors
    for ss, clr, tprior in zip(np.arange(3),
                               ['k', 'b', 'g', 'gray'],
                               [associate_defs.theta_u,
                                associate_defs.theta_c,
                                associate_defs.theta_e,
                                associate_defs.theta_e,
                                ]):
        # Plot
        ax = plt.subplot(gs[ss])
        # Prior
        prior['theta'] = tprior


        if ss < 3:
            # Bayesian + parse
            _, model_mags, model_theta, max_PMix = analysis.run(prior, frb_pre=frb_pre)

            # Distribution
            ptheta = bayesian.pw_Oi(thetas, 1., prior['theta'])
            scl = np.sum(ptheta*dtheta)
            print('scl:', scl)
            convolved_ptheta = None

            # High confidence
            high_conf_theta = []
            nsecure = 0
            for frbA in frb_pre.frbA:
                if np.max(frbA.candidates.P_Ox) > POx_secure:
                    imax = np.argmax(frbA.candidates.P_Ox)
                    high_conf_theta.append(
                        frbA.candidates.iloc[imax].separation /
                        frbA.candidates.iloc[imax].half_light)
                    # Convolve with prior
                    sigR_phi = np.sqrt(frbA.frb.sig_a*frbA.frb.sig_b
                                   ) / frbA.candidates.iloc[imax].half_light
                    conv_theta = convolve(ptheta,
                                          np.exp(-thetas**2/2/sigR_phi**2), mode='full')
                    scl_c = np.sum(conv_theta*dtheta)
                    if convolved_ptheta is None:
                        convolved_ptheta = conv_theta / scl_c
                    else:
                        convolved_ptheta += conv_theta / scl_c
                    # Increment
                    nsecure += 1
            # Nomalize
            high_conf_theta = np.array(high_conf_theta)
            convolved_ptheta /= nsecure

            # Save
            if tprior == associate_defs.theta_e:
                save_high = high_conf_theta.copy()

            # KS test
            #cumsum = np.cumsum(ptheta*dtheta) / scl
            cumsum = np.cumsum(convolved_ptheta*dtheta)
            more_thetas = np.arange(convolved_ptheta.size) * dtheta
            #f_CDF = interp1d(thetas, cumsum)
            f_CDF = interp1d(more_thetas, cumsum)
            res = kstest(high_conf_theta, f_CDF)
            pvalue = res.pvalue

        else:
            # Fit secure with a new exponential
            scale_lengths = np.linspace(0.25, 2., 100)
            P_KS = []
            for scale_length in scale_lengths:
                ptheta = bayesian.pw_Oi(thetas, 1., associate_defs.theta_e,
                                        scale_half=scale_length)
                scl = np.sum(ptheta * dtheta)

                # KS test
                cumsum = np.cumsum(ptheta * dtheta) / scl
                f_CDF = interp1d(thetas, cumsum)
                res = kstest(high_conf_theta, f_CDF)
                # Save
                P_KS.append(res.pvalue)

            imax = np.argmax(P_KS)
            best_sl = scale_lengths[imax]
            print("Best scale length = {}".format(best_sl))
            pvalue = P_KS[imax]
            # Once more
            ptheta = bayesian.pw_Oi(thetas, 1., associate_defs.theta_e, scale_half=best_sl)
            scl = np.sum(ptheta * dtheta)
            prior['theta']['method'] = 'exp '+'{:0.1f}'.format(best_sl)+r'$ \, '+cmhalf+'$'
            #
            high_conf_theta = save_high

        # Plot profiles
        ax.plot(thetas, ptheta/scl,
                label=prior['theta']['method']+' unconvolved', color=clr, alpha=0.3)
        more_thetas = np.arange(convolved_ptheta.size) * dtheta
        ax.plot(more_thetas, convolved_ptheta,
                label=prior['theta']['method']+r' $P_{\rm KS} ='+'{:0.2f}'.format(
                    pvalue)+'$', color=clr)

        # Plot secure
        weights3 = np.ones_like(high_conf_theta) / high_conf_theta.size
        lbl = r'$P(O_i) > '+'{}'.format(POx_secure)+'$ FRBs' if ss == 0 else None
        ax.hist(high_conf_theta, bins=bins_theta, weights=weights3,
                color='darkred', label=lbl,
                histtype='stepfilled')

        # Plot all
        if ss < 3:
            weights2 = np.ones_like(model_theta) / model_theta.size
            lbl = 'All FRB candidates' if ss == 0 else None
            ax.hist(model_theta, weights=weights2, bins=bins_theta,
                    color='darkgrey', label=lbl, histtype='step')



        # Label me
        ax.set_xlabel(r'$\theta/'+cmhalf+'$')
        ax.set_ylabel('PDF')
        ax.xaxis.set_major_locator(plt.MultipleLocator(1.))
        #ax.xaxis.set_major_formatter(FormatStrFormatter(r'$%.3f$'))
        ax.set_xlim(0., 8)
        ax.set_ylim(0., 0.5)

        # Legend
        legend = ax.legend(loc='upper right', scatterpoints=1, borderpad=0.2,
                       handletextpad=handletextpad, fontsize=11.)

        # Font size
        set_fontsize(ax, 13.)


    # End
    plt.tight_layout(pad=0.2, h_pad=0., w_pad=0.1)
    print('Writing {:s}'.format(outfile))
    kwargs = {}
    if 'png' in outfile:
        kwargs['dpi'] = 700
    plt.savefig(outfile, **kwargs)
    plt.close()


def fig_pchance(frb_name='FRB180924', outfile='fig_pchance.png'):
    """

    Args:
        frb_name:
        outfile:
    """
    set_mplrc()
    mpl.rcParams['text.latex.preamble'] = [r'\boldmath']

    if frb_name == 'FRB180924':
        config = frbs.frb180924
    elif frb_name == 'FRB190614':
        config = frbs.frb190614
        config['cut_size'] = 20.
    else:
        raise IOError("Bad FRB")
    config['skip_bayesian'] = True
    frbA = frbassociate.run_individual(config)#, verbose=True, show=True)

    fsz = 9.
    lsz = fsz
    fig = plt.figure(figsize=(5, 5))

    # Image
    x0, y0 = 0.25, 0.1
    aximg = fig.add_axes([x0, y0, 0.95-x0, 0.98-y0], projection=frbA.wcs)

    #blues = plt.get_cmap('Blues')
    greens = plt.get_cmap('Greens')
    vmnx = 200.
    d = aximg.imshow(frbA.hdu.data, cmap=greens, vmin=-vmnx, vmax=vmnx)

    # Label candidates
    for cc, cand in frbA.candidates.iterrows():
        roffset, doffset = 0, 0
        if frb_name == 'FRB180924':
            doffset = 0.9/3600. if cand.separation > 1.5 else -4./3600.
            if cand.separation < 1.5:
                roffset = 10/3600.
            elif np.abs(getattr(cand, frbA.filter)-25.5) < 0.1:
                roffset = 3 / 3600.
            else:
                roffset = 1 / 3600.
        elif frb_name == 'FRB190614':
            doffset = -4./3600. if np.abs(cand.separation-2.18) < 0.1 else 0.9 / 3600.
            roffset = 3/3600 # if np.abs(cand.separation-2.18) < 0.1 else 0 / 3600.

        aximg.text(cand.ra+roffset, cand.dec+doffset,
                   r"$m = {:0.1f}$".format(getattr(cand, frbA.filter))+
                   "\n"+
                   r"$"+cmhalf+"={:0.1f}''$".format(cand.half_light)+
                   "\n"+
                   r"$\theta={:0.1f}''$".format(cand.separation)+
                   "\n"+
                   r"$P^c={:0.2f}$".format(cand.P_c),
                   color='k', transform=aximg.get_transform('icrs'), size=fsz, ha='left')
        #aximg.plot(row['ra'], row['dec'], markers[ss], color='pink',
        #           transform=aximg.get_transform('icrs'), ms=2,
        #           alpha=0.5)

    # FRB localization
    apermap = plot_frb_localization(frbA.frb, frbA.wcs, scale=5.)
    apermap.plot(color='red', lw=1, ls='dotted')

    plt.grid(color='gray', ls='dashed', lw=0.5)
    asz = 11.
    aximg.set_xlabel(r'\textbf{Right Ascension (J2000)}', fontsize=asz)
    aximg.set_ylabel(r'\textbf{Declination (J2000)}', fontsize=asz, labelpad=-1)
    aximg.tick_params(axis='both', labelsize=lsz)

    aximg.set_title(frbA.frb.frb_name)

    # End
    plt.tight_layout(pad=0.2, h_pad=0., w_pad=0.1)
    print('Writing {:s}'.format(outfile))
    kwargs = {}
    if 'png' in outfile:
        kwargs['dpi'] = 700
    plt.savefig(outfile, **kwargs)
    plt.close()


def fig_galaxy_priors(frb_name, outfile='fig_galaxy_priors.png'):
    """

    Args:
        frb_name:
        outfile:

    Returns:

    """
    if frb_name == 'FRB180924':
        config = frbs.frb180924
    elif frb_name == 'FRB190523':
        config = frbs.frb190523
    else:
        raise IOError("Bad FRB")
    config['skip_bayesian'] = True
    frbA = frbassociate.run_individual(config)#, verbose=True, show=True)

    # Pchance
    frbA.calc_pchance()

    # Priors

    # Uniform
    P_U = 0
    frbA.calc_priors(P_U, method='identical')
    frbA.candidates['P_O_u'] = frbA.candidates['P_O'].values.copy()

    # Linear
    #frbA.calc_priors(P_U, method='linear')
    #frbA.candidates['P_O_l'] = frbA.candidates['P_O'].values.copy()

    # Inverse
    frbA.calc_priors(P_U, method='inverse')
    frbA.candidates['P_O_i'] = frbA.candidates['P_O'].values.copy()

    fsz = 9.
    lsz = fsz
    fig = plt.figure(figsize=(3, 6))

    # Image
    x0, y0 = 0.25, 0.5
    aximg = fig.add_axes([x0, y0, 0.95-x0, 0.98-y0], projection=frbA.wcs)

    blues = plt.get_cmap('Blues')
    vmnx = 200.
    d = aximg.imshow(frbA.hdu.data, cmap=blues, vmin=-vmnx, vmax=vmnx)

    # Mark candidates
    markers = ['o', 's', 'D', '^']

    ss = 0
    for _, cand in frbA.candidates.iterrows():
            aximg.plot(cand.ra, cand.dec, markers[ss], color='pink',
                   transform=aximg.get_transform('icrs'), ms=2,
                   alpha=0.5)
            ss += 1

    # FRB localization
    apermap = plot_frb_localization(frbA.frb, frbA.wcs)
    apermap.plot(color='red', lw=1, ls='dotted')

    plt.grid(color='gray', ls='dashed', lw=0.5)
    aximg.set_xlabel('Right Ascension (J2000)', fontsize=fsz)
    aximg.set_ylabel('Declination (J2000)', fontsize=fsz, labelpad=-2)
    aximg.tick_params(axis='both', labelsize=lsz)

    aximg.set_title(frbA.frb.frb_name)

    # ###############################################################
    # Priors
    #
    x0, y0 = 0.18, 0.08
    axp = fig.add_axes([x0, y0, 0.97-x0, 0.47-y0])

    clrs = ['g', 'k', 'b']
    #for lbl, postfix, clr in zip(['Identical', 'Linear', 'Inverse'],
    #                             ['_u', '_l', '_i'],
    #                             clrs):
    for lbl, postfix, clr in zip(['Identical', 'Inverse'],
                                 ['_u', '_i'],
                                 clrs):
        ss = 0
        for _, cand in frbA.candidates.iterrows():
            _lbl = lbl if ss == 0 else None
            axp.scatter(cand.separation, cand['P_O'+postfix],
                    marker=markers[ss], color=clr, label=_lbl)
            ss += 1

    axp.set_xlabel(r'$\theta$ (arcsec)')
    axp.set_ylabel(cPOi)
    axp.set_xlim(0., 10.)
    axp.set_ylim(0., 1.)
    set_fontsize(axp, 9.)

    legend = axp.legend(loc='upper right', scatterpoints=0, borderpad=0.2,
                        handlelength=0, handletextpad=0, fontsize='small')
    # Color the text
    for ss, text in enumerate(legend.get_texts()):
        plt.setp(text, color=clrs[ss])
    # Turn off dots
    for item in legend.legendHandles:
        item.set_visible(False)

    # End
    plt.tight_layout(pad=0.2, h_pad=0., w_pad=0.1)
    print('Writing {:s}'.format(outfile))
    kwargs = {}
    if 'png' in outfile:
        kwargs['dpi'] = 700
    plt.savefig(outfile, **kwargs)
    plt.close()


def fig_pOx_frbs(frb_list=None, outfile='fig_pOx_frbs.png',
                 tmax=None, exp_offset=False,
                 verbose=False):
    """

    Args:
        outfile:

    Returns:

    """
    if frb_list is None:
        frb_list = associate_defs.frb_list
    # Priors
    #priors = [associate_defs.conservative, associate_defs.adopted, associate_defs.extreme]
    priors = [associate_defs.conservative.copy(), associate_defs.adopted.copy()]
    # Adjust priors
    if tmax is not None:
        for prior in priors:
            prior['theta']['max'] = tmax
    if exp_offset:
        for prior in priors:
            prior['theta']['method'] = 'exp'
    #
    clrs = ['g', 'k', 'b']
    Npriors = len(priors)
    width = 0.9/Npriors

    # Plot
    set_mplrc()
    plt.figure(figsize=(15, 5))
    ax = plt.gca()

    ind = np.arange(len(frb_list))
    plt.xticks(ind, frb_list)

    # Loop on FRBs
    for ss, frb_name in enumerate(frb_list):
        config = getattr(frbs, frb_name.lower())
        config['skip_bayesian'] = True
        for pp, prior in enumerate(priors):
            frbA = frbassociate.run_individual(config)
            frbA.calc_priors(prior['U'], method=prior['O'])
            prior['theta']['ang_size'] = frbA.candidates['half_light'].values
            frbA.set_theta_prior(prior['theta'])

            # Calculate p(O_i|x)
            frbA.calc_POx()
            # Reverse sort
            frbA.candidates = frbA.candidates.sort_values('P_Ox', ascending=False)

            cum_POx = 0
            for cc, cand in frbA.candidates.iterrows():
                # Bar me!
                ax.bar(ss-0.45 + (pp+0.5)*width, cand.P_Ox, width, color=clrs[pp],
                       bottom=cum_POx, edgecolor=clrs[pp], fill=(cum_POx==0))
                cum_POx += cand.P_Ox
            # Print
            if verbose:
                print(frb_name+'\n', frbA.candidates[['id', frbA.filter, 'half_light',
                                  'separation', 'P_O', 'P_Ox']])

    # Label me
    #ax.set_xlabel(r'$\theta/\theta_{1/2}$')
    ax.set_ylabel(cPOxi)
    #ax.xaxis.set_major_locator(plt.MultipleLocator(1.))
    #ax.xaxis.set_major_formatter(FormatStrFormatter(r'$%.3f$'))
    ax.set_ylim(0., 1.)

    # Legend
    legend = ax.legend(loc='upper right', scatterpoints=1, borderpad=0.2,
                       handletextpad=handletextpad, fontsize='large')

    set_fontsize(ax, 13.)
    # End
    plt.tight_layout(pad=0.2, h_pad=0., w_pad=0.1)
    print('Writing {:s}'.format(outfile))
    kwargs = {}
    if 'png' in outfile:
        kwargs['dpi'] = 700
    plt.savefig(outfile, **kwargs)
    plt.close()


def fig_vary_PU(frb_list=None, outfile='fig_vary_PU.png'):
    """

    Args:
        outfile:

    Returns:

    """
    sns.set_theme()
    sns.set_style('whitegrid')
    sns.set_context('paper')

    if frb_list is None:
        frb_list = associate_defs.frb_list
    # Priors
    aprior = associate_defs.adopted.copy()

    # Loop on PU
    #PUs = [0., 0.1, 0.5, 0.9]
    PUs = np.linspace(0., 0.1, 10)

    POx_frbs = np.zeros((len(frb_list), len(PUs)))
    for ss, frb_name in enumerate(frb_list):
        config = getattr(frbs, frb_name.lower())
        config['skip_bayesian'] = True
        for pp, PU in enumerate(PUs):
            frbA = frbassociate.run_individual(config)
            frbA.calc_priors(PU, method=aprior['O'])
            aprior['theta']['ang_size'] = frbA.candidates['half_light'].values
            frbA.set_theta_prior(aprior['theta'])

            # Calculate p(O_i|x)
            frbA.calc_POx()
            # Reverse sort
            frbA.candidates = frbA.candidates.sort_values('P_Ox', ascending=False)

            # Save
            POx_frbs[ss, pp] = frbA.candidates.P_Ox.values[0]
            #if PU > 0.:
                #embed(header='613 of figs')
                #frbA.candidates[['ra', 'dec', 'P_O', 'P_Ox', 'p_xO']]

    # Plot
    set_mplrc()
    plt.figure(figsize=(15, 5))
    ax = plt.gca()

    # Colors
    cm = plt.get_cmap('jet')
    N = len(frb_list) + 1
    plt.rcParams["axes.prop_cycle"] = plt.cycler("color", cm(np.linspace(0, 1, N)))

    # Plot
    for kk, frb_name in enumerate(frb_list):
        sns.lineplot(x=PUs, y=POx_frbs[kk,:], ax=ax, label=frb_name)

    # Label me
    ax.set_xlabel(r'$P(U)$')
    ax.set_ylabel(cPOxi)
    #ax.xaxis.set_major_locator(plt.MultipleLocator(1.))
    #ax.xaxis.set_major_formatter(FormatStrFormatter(r'$%.3f$'))
    ax.set_ylim(0., 1.)

    # Legend
    legend = ax.legend(loc='upper right', scatterpoints=1, borderpad=0.2,
                       handletextpad=handletextpad, fontsize='large')

    set_fontsize(ax, 13.)
    # End
    plt.tight_layout(pad=0.2, h_pad=0., w_pad=0.1)
    print('Writing {:s}'.format(outfile))
    kwargs = {}
    if 'png' in outfile:
        kwargs['dpi'] = 700
    plt.savefig(outfile, **kwargs)
    plt.close()



def fig_frb_excess(frb_list=None, outfile='fig_frb_excess.png', verbose=False):
    """

    Args:
        outfile:

    Returns:

    """
    sns.set_theme()
    sns.set_style('whitegrid')
    sns.set_context('paper')

    if frb_list is None:
        frb_list = associate_defs.frb_list

    # Plot
    plt.figure(figsize=(8, 8))
    ax = plt.gca()

    # Loop on FRBs
    all_seps = []
    for ss, frb_name in enumerate(frb_list):
        config = getattr(frbs, frb_name.lower())
        config['skip_bayesian'] = True
        # Go to 30''
        max = 60.
        config['max_radius'] = max
        config['cut_size'] = max
        config['cand_separation'] = max * units.arcsec
        #
        try:
            frbA = frbassociate.run_individual(config)
        except:
            embed(header='fail')
        # Save
        all_seps += frbA.candidates.separation.values.tolist()

    # Histogram
    _ = sns.histplot(x=all_seps, ax=ax)#, log_scale=(False,True))

    ax.set_xlabel(r'Separation (arcsec)')
    ax.set_ylabel(r'Number')

    # Font size
    set_fontsize(ax, 15.)

    # End
    print('Writing {:s}'.format(outfile))
    kwargs = {}
    if 'png' in outfile:
        kwargs['dpi'] = 700
    plt.savefig(outfile, **kwargs)
    plt.close()


def fig_mag_vs_DM(outfile='fig_mag_vs_DM.png'):
    """
    mini Maquart relation

    Args:
        outfile:

    Returns:

    """
    set_mplrc()

    # Load f_mL
    f_mL = analysis.load_f_mL()

    # Init prior
    prior = associate_defs.adopted.copy()

    # Plot
    plt.figure(figsize=(8, 5))
    gs = gridspec.GridSpec(1,1)

    cm = plt.get_cmap('jet')

    # Bayesian + parse
    frbA_tbl, model_mags, model_theta, max_PMix = analysis.run(prior)

    # Colors
    N = len(frbA_tbl) + 1
    plt.rcParams["axes.prop_cycle"] = plt.cycler("color", cm(np.linspace(0, 1, N)))

    # Plot -- must come after colors
    ax = plt.subplot(gs[0])

    # Plot em
    Pmin = 0.01
    unit_size = 150.
    order = [0,1,6,3,2,4,5,7,8,9,10,11,12]
    DMcosmic = [frb_row.frbA.frb.DM.value - frb_row.frbA.frb.DMISM.value - 100
                for _, frb_row in frbA_tbl.iterrows()]
    order = np.argsort(DMcosmic)
    #for _, frb_row in frbA_tbl.iterrows():
    for ss in order:
        # Restrict to candiates with min
        frb_row = frbA_tbl.iloc[ss]
        ok_c = frb_row.frbA.candidates.P_Ox > Pmin
        N_c = np.sum(ok_c)
        # Scatter plot em
        ax.scatter([frb_row.frbA.frb.DM.value-frb_row.frbA.frb.DMISM.value-100]*N_c,
                   frb_row.frbA.candidates[ok_c][frb_row.frbA.filter],
                   label=frb_row.frb,
                   edgecolor='darkslategrey',
                   marker='o',
                   s=unit_size * frb_row.frbA.candidates[ok_c].P_Ox,
                   )

    # Fiducial relation
    L_Lstar = 0.270           # Grabbed from PL
    numL = int(np.round(1./L_Lstar))
    log10_L_Lstar = np.log10(L_Lstar)
    std_log10_L_Lstar = 0.47  # Grabbed from fig_PL

    # Load up DM_cosmic
    dm_cosmic, z = igm.average_DM(1., cumul=True)
    f_DMz = interp1d(z, dm_cosmic.value)

    zval = np.linspace(0.02, 1., 100)
    DMs = f_DMz(zval)

    mLstar = f_mL(zval)
    m_FRB = mLstar - 2.5*log10_L_Lstar

    ax.plot(DMs, m_FRB, 'k-', label=r'$L=L*/${}'.format(numL))

    # Label me
    ax.set_xlabel(r'DM$_{\rm cosmic} \equiv$ DM$_{\rm FRB}$ - DM$_{\rm ISM}$ - 100')
    ax.set_ylabel(r'$m_r$')
    #ax.xaxis.set_major_locator(plt.MultipleLocator(1.))
    #ax.xaxis.set_major_formatter(FormatStrFormatter(r'$%.3f$'))
    ax.set_xlim(0., 800.)

    # Legend
    legend = ax.legend(loc='lower right', scatterpoints=1, borderpad=0.0,
                   handletextpad=handletextpad, fontsize=10.)

    # Font size
    set_fontsize(ax, 15.)

    # End
    plt.tight_layout(pad=0.2, h_pad=0., w_pad=0.1)
    print('Writing {:s}'.format(outfile))
    kwargs = {}
    if 'png' in outfile:
        kwargs['dpi'] = 700
    plt.savefig(outfile, **kwargs)
    plt.close()


def fig_PL(outfile='fig_PL.png'):
    """

    Args:
        outfile:

    Returns:

    """
    set_mplrc()

    sns.set_theme()
    sns.set_style('whitegrid')
    sns.set_context('paper')

    # Load f_mL
    f_mL = analysis.load_f_mL()

    frb_pre = analysis.get_candidates()

    # Init prior
    prior = associate_defs.adopted.copy()

    # Plot
    plt.figure(figsize=(4, 8))
    gs = gridspec.GridSpec(2,1)


    # Bayesian + parse
    _, model_mags, model_theta, max_PMix = analysis.run(prior, frb_pre=frb_pre)

    #bins_L = np.linspace(0.01, 10., 20)
    bins_L = np.linspace(-2., 1, 15)
    #bins_L = np.linspace(0.1, 2., 20)
    #bins_L = np.linspace(16., 25., 20) # Testing only

    # High confidence
    high_conf_L, high_conf_z, high_conf_mr = [], [], []
    for frbA in frb_pre.frbA:
        if np.max(frbA.candidates.P_Ox) > associate_defs.POx_secure:
            # This is klunky
            hg = frbA.frb.grab_host()
            if 'z_spec' not in hg.redshift.keys():
                continue
            imax = np.argmax(frbA.candidates.P_Ox)
            # Magnitude
            m = frbA.candidates.iloc[imax][frbA.filter]
            # m_r(L*)
            m_r_Lstar = float(f_mL(frbA.frb.z))
            # Now magnitude time
            log10_Lstar = (m_r_Lstar-m)/2.5
            print('FRB, z, m, m_L*, log10_L: ', frbA.frb.frb_name, frbA.frb.z, m,
                  m_r_Lstar, log10_Lstar)
            # Save
            high_conf_L.append(log10_Lstar)
            high_conf_z.append(frbA.frb.z)
            high_conf_mr.append(m)
    high_conf_L = np.array(high_conf_L)
    high_conf_z = np.array(high_conf_z)
    high_conf_mr = np.array(high_conf_mr)

    # Plot m_r vs. z
    ax = plt.subplot(gs[0])
    sns.scatterplot(x=high_conf_z, y=high_conf_mr, ax=ax)
    ax.set_xlabel(r'$z$')
    ax.set_ylabel(r'$m_r$')

    # Add a line
    zs = np.linspace(0.02, 0.5, 100)
    mrss = f_mL(zs)
    ax.plot(zs, mrss, 'k--')

    fsz= 15.
    set_fontsize(ax, fsz)

    # Plot
    ax = plt.subplot(gs[1])

    weights = np.ones_like(high_conf_L) / high_conf_L.size
    lbl = r'$P(O_i) > '+'{}'.format(associate_defs.POx_secure)+'$ FRBs'
    ax.hist(high_conf_L, bins=bins_L, weights=weights,
            color='b', label=lbl,
            histtype='stepfilled')

    # Stats
    print("Median L/L* = {}".format(np.median(10**high_conf_L)))
    print("RMS log(L/L*) = {}".format(np.std(high_conf_L)))

    # Label me
    ax.set_xlabel(r'$\log_{10} \, (L/L*)$')
    ax.set_ylabel('PDF')
    ax.xaxis.set_major_locator(plt.MultipleLocator(1.))
    #ax.xaxis.set_major_formatter(FormatStrFormatter(r'$%.3f$'))
    ax.set_ylim(0., 0.3)

    # Legend
    legend = ax.legend(loc='upper right', scatterpoints=1, borderpad=0.2,
                   handletextpad=handletextpad, fontsize=13.)

    # Font size
    set_fontsize(ax, fsz)

    # End
    plt.tight_layout(pad=0.2, h_pad=0., w_pad=0.1)
    print('Writing {:s}'.format(outfile))
    kwargs = {}
    if 'png' in outfile:
        kwargs['dpi'] = 700
    plt.savefig(outfile, **kwargs)
    plt.close()



def fig_Pcolor(outfile='fig_Pcolor.png'):
    """

    Args:
        outfile:

    Returns:

    """
    set_mplrc()

    # Load f_mL
    f_mL = analysis.load_f_mL()

    frb_pre = analysis.get_candidates()

    # Init prior
    prior = associate_defs.adopted.copy()

    # Plot
    plt.figure(figsize=(7, 7))


    # Bayesian + parse
    _, model_mags, model_theta, max_PMix = analysis.run(prior, frb_pre=frb_pre)

    #bins_L = np.linspace(0.01, 10., 20)
    bins_L = np.linspace(-2., 1, 15)
    #bins_L = np.linspace(0.1, 2., 20)
    #bins_L = np.linspace(16., 25., 20) # Testing only

    # High confidence
    high_conf_Mstar, high_conf_ur = [], []
    for frbA in frb_pre.frbA:
        if np.max(frbA.candidates.P_Ox) > associate_defs.POx_secure:
            # This is klunky
            hg = frbA.frb.grab_host()
            if 'z_spec' not in hg.redshift.keys():
                continue
            imax = np.argmax(frbA.candidates.P_Ox)
            # u-r
            if 'u-r' not in hg.derived.keys():
                continue
            ur = hg.derived['u-r']
            Mstar = hg.derived['Mstar']
            # Save
            high_conf_Mstar.append(np.log10(Mstar))
            high_conf_ur.append(ur)
    high_conf_Mstar = np.array(high_conf_Mstar)
    high_conf_ur = np.array(high_conf_ur)

    # Plot m_r vs. z
    df = pandas.DataFrame(dict(Mstar=high_conf_Mstar, ur=high_conf_ur))
    jg = sns.jointplot(data=df, x='Mstar', y='ur')
    jg.ax_marg_x.set_xlim(8, 10.5)
    jg.ax_marg_y.set_ylim(0.5, 2.0)
    jg.ax_joint.set_xlabel(r'$\log_{10} \, (M*/M_\odot)$')
    jg.ax_joint.set_ylabel(r'$u-r$')
    jg.ax_joint.minorticks_on()

    jg.ax_joint.yaxis.set_major_locator(plt.MultipleLocator(0.5))
    jg.ax_joint.xaxis.set_major_locator(plt.MultipleLocator(1.0))
    #jg.ax_joint.yaxis.set_major_formatter(FormatStrFormatter(r'$%.3f$'))

    # Font size
    set_fontsize(jg.ax_joint, 15.)

    # End
    plt.tight_layout(pad=0.2, h_pad=0., w_pad=0.1)
    print('Writing {:s}'.format(outfile))
    kwargs = {}
    if 'png' in outfile:
        kwargs['dpi'] = 700
    plt.savefig(outfile, **kwargs)
    plt.close()


def fig_Phalflight(outfile='fig_Phalflight.png'):
    """

    Args:
        outfile:

    Returns:

    """
    set_mplrc()

    # Load f_mL
    f_mL = analysis.load_f_mL()

    frb_pre = analysis.get_candidates()

    # Init prior
    prior = associate_defs.adopted.copy()

    # Plot
    plt.figure(figsize=(7, 7))


    # Bayesian + parse
    _, model_mags, model_theta, max_PMix = analysis.run(prior, frb_pre=frb_pre)

    #bins_L = np.linspace(0.01, 10., 20)
    bins_L = np.linspace(-2., 1, 15)
    #bins_L = np.linspace(0.1, 2., 20)
    #bins_L = np.linspace(16., 25., 20) # Testing only

    # High confidence
    high_conf_z, high_conf_hl = [], []
    for frbA in frb_pre.frbA:
        if np.max(frbA.candidates.P_Ox) > associate_defs.POx_secure:
            # Half light
            imax = np.argmax(frbA.candidates.P_Ox)
            # Save
            high_conf_z.append(frbA.frb.z)
            high_conf_hl.append(frbA.candidates.iloc[imax].half_light)
    high_conf_z = np.array(high_conf_z)
    high_conf_hl = np.array(high_conf_hl)

    # Plot m_r vs. z
    df = pandas.DataFrame(dict(z=high_conf_z, hl=high_conf_hl))
    jg = sns.jointplot(data=df, x='z', y='hl')
    jg.ax_marg_x.set_xlim(0, 0.6)
    #jg.ax_marg_y.set_ylim(0.5, 2.5)
    jg.ax_joint.set_xlabel(r'$z$')
    jg.ax_joint.set_ylabel(chalf+' (arcsec)')

    #jg.ax_joint.yaxis.set_major_locator(plt.MultipleLocator(0.5))
    #jg.ax_joint.yaxis.set_major_formatter(FormatStrFormatter(r'$%.3f$'))

    # Font size
    #set_fontsize(ax, 13.)

    # End
    plt.tight_layout(pad=0.2, h_pad=0., w_pad=0.1)
    print('Writing {:s}'.format(outfile))
    kwargs = {}
    if 'png' in outfile:
        kwargs['dpi'] = 700
    plt.savefig(outfile, **kwargs)
    plt.close()


def fig_Psizesep(outfile='fig_Psizesep.png'):
    """

    Args:
        outfile:

    Returns:

    """
    set_mplrc()

    sns.set_theme()
    sns.set_style('whitegrid')
    sns.set_context('paper')

    # Load f_mL
    f_mL = analysis.load_f_mL()

    frb_pre = analysis.get_candidates()

    # Init prior
    prior = associate_defs.adopted.copy()

    # Plot
    plt.figure(figsize=(7, 7))


    # Bayesian + parse
    _, model_mags, model_theta, max_PMix = analysis.run(prior, frb_pre=frb_pre)

    #bins_L = np.linspace(0.01, 10., 20)
    bins_L = np.linspace(-2., 1, 15)
    #bins_L = np.linspace(0.1, 2., 20)
    #bins_L = np.linspace(16., 25., 20) # Testing only

    # High confidence
    high_conf_size, high_conf_sep = [], []
    max_sep = 0.
    for frbA in frb_pre.frbA:
        if np.max(frbA.candidates.P_Ox) > associate_defs.POx_secure:
            # Max
            imax = np.argmax(frbA.candidates.P_Ox)
            # Measures
            ang_phys = associate_defs.cosmo.kpc_proper_per_arcmin(frbA.frb.z)
            high_conf_sep.append((frbA.candidates.iloc[imax].separation*units.arcsec *
                                 ang_phys).to('kpc').value)
            high_conf_size.append(
                ((frbA.candidates.iloc[imax].half_light*units.arcsec) *
                 ang_phys).to('kpc').value)
            #
            if high_conf_sep[-1] > max_sep:
                max_sep = high_conf_sep[-1]
                max_FRB = frbA.frb
    high_conf_size = np.array(high_conf_size)
    high_conf_sep = np.array(high_conf_sep)
    print("Max sep: {}".format(max_FRB))

    # Plot m_r vs. z
    df = pandas.DataFrame(dict(size=high_conf_size, sep=high_conf_sep))
    jg = sns.jointplot(data=df, x='sep', y='size')
    #jg.ax_marg_x.set_xlim(0, 0.6)
    jg.ax_marg_y.set_ylim(0., 7)
    jg.ax_joint.set_xlabel('Physical Separation (kpc)')
    jg.ax_joint.set_ylabel('Galaxy size (kpc)')

    # One-to-one line
    jg.ax_joint.plot([0., 11.], [0., 11.], 'k--')

    #jg.ax_joint.yaxis.set_major_locator(plt.MultipleLocator(0.5))
    #jg.ax_joint.yaxis.set_major_formatter(FormatStrFormatter(r'$%.3f$'))

    # Font size
    set_fontsize(jg.ax_joint, 15.)

    # End
    plt.tight_layout(pad=0.2, h_pad=0., w_pad=0.1)
    print('Writing {:s}'.format(outfile))
    kwargs = {}
    if 'png' in outfile:
        kwargs['dpi'] = 700
    plt.savefig(outfile, **kwargs)
    plt.close()


def fig_segm(outfile='fig_segm.png'):
    """

    Args:
        outfile:

    Returns:

    """
    set_mplrc()

    # Plot
    plt.figure(figsize=(9, 5))
    gs = gridspec.GridSpec(1,2)

    # FRB 180924
    ax0 = plt.subplot(gs[0])
    config0 = frbs.frb180924
    frbA_180924 = frbassociate.run_individual(config0)

    cmap = frbA_180924.segm.make_cmap()
    ax0.imshow(frbA_180924.segm, origin='lower', cmap=cmap, interpolation='nearest',
               extent=[-config0['cut_size']/2, config0['cut_size'] / 2,
                       -config0['cut_size'] / 2, config0['cut_size'] / 2])
    ax0.set_title('FRB 180924')

    # FRB localization
    apermap = plot_frb_localization(frbA_180924.frb,
                                    frbA_180924.wcs)
    apermap.plot(color='white', lw=1)

    # FRB 190523
    ax1 = plt.subplot(gs[1])
    config1 = frbs.frb190523
    frbA_190523 = frbassociate.run_individual(config1)

    cmap = frbA_190523.segm.make_cmap()
    ax1.imshow(frbA_190523.segm, origin='lower', cmap=cmap, interpolation='nearest',
               extent = [-config1['cut_size'] / 2, config1['cut_size'] / 2,
              -config1['cut_size'] / 2, config1['cut_size'] / 2])
    ax1.set_title('FRB 190523')

    apermap = plot_frb_localization(frbA_190523.frb,
                                    frbA_190523.wcs)
    apermap.plot(color='white', lw=1)

    # Label me
    for ax in [ax0, ax1]:
        ax.set_xlabel('arcsec')
        ax.set_ylabel('arcsec')
    #ax.xaxis.set_major_locator(plt.MultipleLocator(1.))
    #ax.xaxis.set_major_formatter(FormatStrFormatter(r'$%.3f$'))

    # Font size
    for ax in [ax0, ax1]:
        set_fontsize(ax, 15.)

    # End
    plt.tight_layout(pad=0.2, h_pad=0., w_pad=0.1)
    print('Writing {:s}'.format(outfile))
    kwargs = {}
    if 'png' in outfile:
        kwargs['dpi'] = 700
    plt.savefig(outfile, **kwargs)
    plt.close()


def fig_galaxy_posteriors(frb_name, outfile='fig_galaxy_posteriors.png',
                          raoffsets=None, decoffsets=None, cmap=None,
                          frb_color='r', vmnx=(None,None), cut_low=False,
                          add_faint=None, prior_mode='adopted',
                          min_POx=0.01):
    """

    Args:
        frb_name:
        outfile:
        add_faint (tuple, optional):
            separation, half-light, mag

    Returns:

    """
    if frb_name in ['FRB180301']:
        config = analysis.frb180301
    else:
        config = getattr(frbs, frb_name.lower())
    # Run
    config['skip_bayesian'] = True
    frbA = frbassociate.run_individual(config)#, verbose=True, show=True)

    # Add fake?
    if add_faint:
        fake_coord = frbA.frb.coord.directional_offset_by(0. * units.deg,
                                                          add_faint[0] * units.arcsec)
        fake_row = frbA.candidates.iloc[0]
        # Fill in
        fake_row.coords = fake_coord
        fake_row.ra = fake_coord.ra.value
        fake_row.dec = fake_coord.dec.value
        fake_row.separation = add_faint[0]
        fake_row.half_light = add_faint[1]
        fake_row[frbA.filter] = add_faint[2]
        fake_row.id = 99
        # Add it in
        frbA.candidates = frbA.candidates.append(fake_row, ignore_index=True)
        # Redo chance
        frbA.calc_pchance()

    # Prior
    prior = getattr(associate_defs, prior_mode).copy()

    # Set priors
    frbA.calc_priors(prior['U'], method=prior['O'])
    prior['theta']['ang_size'] = frbA.candidates['half_light'].values
    frbA.set_theta_prior(prior['theta'])

    # Calculate p(O_i|x)
    frbA.calc_POx()
    frbA.candidates = frbA.candidates.sort_values('P_Ox', ascending=False)

    # Plot
    fsz = 10.
    lsz = fsz
    fig = plt.figure(figsize=(5, 5))

    # Image
    x0, y0 = 0.25, 0.1
    aximg = fig.add_axes([x0, y0, 0.95-x0, 0.98-y0], projection=frbA.wcs)

    if cmap is None:
        cmap = plt.get_cmap('Blues')
    d = aximg.imshow(frbA.hdu.data, cmap=cmap, vmin=vmnx[0], vmax=vmnx[1])

    # Mark candidates
    markers = ['o', 's', 'D', '^']


    # Label candidates
    ss = 0
    for cc, cand in frbA.candidates.iterrows():
        if cut_low and cand.P_Ox < min_POx:
            continue
        # Offset
        roffset = 0. if raoffsets is None else raoffsets[ss]/3600.
        doffset = 0. if decoffsets is None else decoffsets[ss]/3600.
        #
        aximg.text(cand.ra+roffset, cand.dec+doffset,
                   r"$P(O)"+"={:0.2f}$".format(cand.P_O)+
                   "\n"+
                   r"$P(O_i|x) = {:0.2f}$".format(cand.P_Ox),
                   color='k', transform=aximg.get_transform('icrs'), size=fsz, ha='left')
        #
        ss += 1


    # FRB localization
    apermap = plot_frb_localization(frbA.frb, frbA.wcs)
    apermap.plot(color=frb_color, lw=1, ls='dashed')

    plt.grid(color='gray', ls='dashed', lw=0.5)
    aximg.set_xlabel('Right Ascension (J2000)', fontsize=fsz)
    aximg.set_ylabel('Declination (J2000)', fontsize=fsz, labelpad=-1)
    aximg.tick_params(axis='both', labelsize=lsz)

    aximg.set_title(frbA.frb.frb_name)

    # End
    #plt.tight_layout(pad=0.2, h_pad=0., w_pad=0.1)
    print('Writing {:s}'.format(outfile))
    kwargs = {}
    if 'png' in outfile:
        kwargs['dpi'] = 700
    plt.savefig(outfile, **kwargs)
    plt.close()



def plot_frb_localization(ifrb, wcs, scale=1.):
    # FRB localization
    theta = ifrb.eellipse['theta']
    aper = SkyEllipticalAperture(positions=ifrb.coord,
                                 a=scale*ifrb.sig_a * units.arcsecond,
                                 b=scale*ifrb.sig_b * units.arcsecond,
                                 theta=theta * units.deg)
    apermap = aper.to_pixel(wcs)
    return apermap


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
        fig_theta_priors()

    # Priors
    if flg_fig & (2**1):
        fig_galaxy_priors('FRB180924', 'fig_priors_180924.png')
        fig_galaxy_priors('FRB190523', 'fig_priors_190523.png')

    # P(O|x)
    if flg_fig & (2**2):
        fig_pOx_frbs()
        #fig_pOx_frbs(tmax=10., outfile='fig_pOx_frbs_tmax10.png')
        #fig_pOx_frbs(exp_offset=True, outfile='fig_pOx_frbs_exp.png')

    # PChance
    if flg_fig & (2**3):
        fig_pchance()

    # theta posterior vs. prior
    if flg_fig & (2**4):
        fig_prior_vs_posterior()

    # mag vs. DM
    if flg_fig & (2**5):
        fig_mag_vs_DM()

    # p(L)
    if flg_fig & (2**6):
        fig_PL()

    # p(L)
    if flg_fig & (2**7):
        fig_segm()

    # Posteriors
    if flg_fig & (2**8):
        fig_galaxy_posteriors('FRB180924', 'fig_posteriors_180924_conservative.png',
                              raoffsets=(5., 2.5, 2., 3.), prior_mode='conservative',
                              decoffsets=(-5.5, 1.5, 0.9, -3), vmnx=(-200,200.),
                              cmap='Reds')
        fig_galaxy_posteriors('FRB180924', 'fig_posteriors_180924.png',
                              cmap='Greens',
                              raoffsets=(5., 2.5, 2., 3.),
                              decoffsets=(-5.5, 1.5, 0.9, -3), vmnx=(-200,200.))
        fig_galaxy_posteriors('FRB121102', 'fig_posteriors_121102.png',
                              vmnx=(6075, 6426),
                              cut_low=True)
        fig_galaxy_posteriors('FRB190614', 'fig_posteriors_190614.png',
                              cmap='Reds', frb_color='k',
                              vmnx=(-247, 345),
                              decoffsets=(1.0, -3),
                              cut_low=True)
        fig_galaxy_posteriors('FRB190523', 'fig_posteriors_190523.png',
                              cmap='Reds', frb_color='k',
                              raoffsets=(8., 20.5, 7., 3.),
                              decoffsets=(2., -3.5, 2., -3), vmnx=(-200,200.))
        fig_galaxy_posteriors('FRB181112', 'fig_posteriors_181112.png',
                              cmap='Reds', frb_color='k', vmnx=(7700,9500),
                              cut_low=True,
                              raoffsets=(0., -3., 7., 3.),
                              decoffsets=(-3., -1., 2., -3))
        fig_galaxy_posteriors('FRB190611', 'fig_posteriors_190611.png',
                              cmap='Reds', frb_color='k', vmnx=(4833,5146),
                              cut_low=True, add_faint=(1., 0.3, 25.5),
                              #raoffsets=(0., -3., 7., 3.),
                              #decoffsets=(-3., -1., 2., -3))
                              )

    # P(color)
    if flg_fig & (2**9):
        fig_Pcolor()

    # P(half light)
    if flg_fig & (2**10):
        fig_Phalflight()

    # FRB excess
    if flg_fig & (2**11):
        fig_frb_excess()

    # P(size,sep)
    if flg_fig & (2**12):
        fig_Psizesep()

    # Vary P(U)
    if flg_fig & (2**13):
        fig_vary_PU()


# Command line execution
if __name__ == '__main__':


    if len(sys.argv) == 1:
        flg_fig = 0
        #flg_fig += 2**0   # P(w|O)
        #flg_fig += 2**1   # P(O)
        #flg_fig += 2**2   # p(O|x) for FRBs
        flg_fig += 2**3   # PChance
        #flg_fig += 2**4   # Prior vs. posterior for p(w|O)
        #flg_fig += 2**5   # m-Macquart mag vs. DM
        #flg_fig += 2**6   # P(L)
        #flg_fig += 2**7   # Segmentation images
        #flg_fig += 2**8   # Individual FRB Posterior
        #flg_fig += 2**9   # P(color)
        #flg_fig += 2**10   # P(half_light)
        #flg_fig += 2**11   # FRB excess
        #flg_fig += 2**12   # P(sizes, sep)
        #flg_fig += 2**13   # Vary P(U)
    else:
        flg_fig = sys.argv[1]

    main(flg_fig)