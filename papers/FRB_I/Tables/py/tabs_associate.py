#Module for Tables for the Bhandari+19 Hosts paper
# Imports
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os, sys
import copy

import pandas

from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units

from frb.associate import frbassociate
from frb.associate import frbs
from frb import frb

# Local
sys.path.append(os.path.abspath("../Analysis/py"))
import associate_defs

from IPython import embed

# Summary table of results
def mktab_priors(outfile='tab_priors.tex', sub=False):

    priors = [associate_defs.conservative, associate_defs.adopted]#, associate_defs.extreme]

    # Open
    tbfil = open(outfile, 'w')

    # Header
    #tbfil.write('\\clearpage\n')
    tbfil.write('\\begin{deluxetable}{cccccccccccccccc}\n')
    #tbfil.write('\\rotate\n')
    tbfil.write('\\tablewidth{0pc}\n')
    tbfil.write('\\tablecaption{Prior Sets\\label{tab:priors}}\n')
    tbfil.write('\\tabletypesize{\\footnotesize}\n')
    tbfil.write('\\tablehead{\\colhead{Set}  \n')
    tbfil.write('& \\colhead{\\PO} & \\colhead{\\PU} & \\colhead{\\poffset}\n')
    tbfil.write('& \\colhead{\\thmax/\\halflight}\n')
    #tbfil.write("& (deg) & (deg) & ($''$) & & (deg) & (deg) & ($''$) & ($''$) & ($''$) & ($''$) & (mag)\n")
    #tbfil.write("\\\\ (1) & (2) & (3) & (4) & (5) & (6) & (7) & (8) & (9) & (10) & (11) & (12) & (13) & (14) & (15)")
    tbfil.write('} \n')

    tbfil.write('\\startdata \n')

    # Loop on priors
    for pp, prior in enumerate(priors):
        sline = ''

        # Name
        sline += prior['name']

        # P(O)
        sline += '&' + prior['O']

        # P(U)
        if prior['U'] >= 0:
            if prior['U'] == 0:
                sline += '& {}'.format(int(prior['U']))
            else:
                sline += '& {:0.2f}'.format(prior['U'])
        else:
            sline += '&' + 'chance'

        # P(theta)
        sline += '&' + prior['theta']['method']

        # theta_max
        sline += '& {}'.format(int(prior['theta']['max']))

        tbfil.write(sline + '\\\\ \n')

    # End end
    tbfil.write('\\hline \n')


    tbfil.write('\\enddata \n')

    #tbfil.write('\\tablenotetext{a}{Spectroscopic redshifts are reported to 4 significant digits.  Photometric to 2.} \n')
    #tbfil.write('the gas below the line $\\rm DEC_{\\rm off} = \\aslope RA_{\\rm off} \\ayoff$}\n')
    #tbfil.write('\\tablecomments{Column 1: FRB source. Columns 2 and 3: R.A. and Decl. of the FRB (J2000). Column 4: FRB error ellipse. Column 5: FRB classication. Repeating = yes(y)/no(n). Column 6 and 7: R.A. and Dec. of the associated host galaxy (J2000). Column 8: projected angular offset of the FRB to the host galaxy center. Column 9: association radius $\delta x$ \citep{Tunnicliffe14}. Column 10: angular effective radius of the host measured from a sersic model using GALFIT \citep{galfit} on the $i$-band images (or equivalent). Column 11: effective search radius \citep{Bloom02}. Column 12: measured apparent magnitude of the host. Column 13: filter used for the magnitude measurement. Column 14: probability of chance coincidence using the \citet{Bloom02} formalism. Column 15: sample designations following the criteria outlined in $\S$~\\ref{ssec:associate}.}\n')
    # End
    tbfil.write('\\end{deluxetable} \n')

    tbfil.close()
    print('Wrote {:s}'.format(outfile))


# Summary table of results
def mktab_frbs(outfile='tab_frbs.tex', sub=False):


    # Open
    tbfil = open(outfile, 'w')

    # Header
    #tbfil.write('\\clearpage\n')
    tbfil.write('\\begin{deluxetable*}{cccccccccccccccc}\n')
    #tbfil.write('\\rotate\n')
    tbfil.write('\\tablewidth{0pc}\n')
    tbfil.write('\\tablecaption{FRBs Analyzed\\label{tab:frbs}}\n')
    tbfil.write('\\tabletypesize{\\footnotesize}\n')
    tbfil.write('\\tablehead{\\colhead{FRB}  \n')
    tbfil.write('& \\colhead{\\rafrb} & \\colhead{\\decfrb}\n')
    tbfil.write('& \\colhead{\\eea} & \\colhead{\\eeb} & \\colhead{\\eep}\n')
    tbfil.write('& \\colhead{Filter}\n')
    tbfil.write('\\\\')
    tbfil.write("& (deg) & (deg) & ($''$) & ($''$) & (deg) \n")
    #tbfil.write("\\\\ (1) & (2) & (3) & (4) & (5) & (6) & (7) & (8) & (9) & (10) & (11) & (12) & (13) & (14) & (15)")
    tbfil.write('} \n')

    tbfil.write('\\startdata \n')

    # Loop on priors
    for pp, frb_name in enumerate(associate_defs.frb_list):
        # Load FRB
        FRB = frb.FRB.by_name(frb_name)
        config = getattr(frbs, frb_name.lower())

        sline = ''

        # Name
        sline += frb_name.upper()

        # RA
        sline += '& {:0.5f}'.format(FRB.coord.ra.value)

        # Dec
        sline += '& {:0.5f}'.format(FRB.coord.dec.value)

        # Error ellipse a  (arcsec)
        sline += '& {:0.2f}'.format(FRB.sig_a)

        # Error ellipse b  (arcsec)
        sline += '& {:0.2f}'.format(FRB.sig_b)

        # Error ellipse PA  (deg)
        sline += '& {:0.1f}'.format(FRB.eellipse['theta'])

        # Filter
        sline += '& {}'.format(config['filter'].replace('_', '\_'))

        tbfil.write(sline + ' \\\\ \n')

    # End end
    tbfil.write('\\hline \n')


    tbfil.write('\\enddata \n')

    #tbfil.write('\\tablenotetext{a}{Spectroscopic redshifts are reported to 4 significant digits.  Photometric to 2.} \n')
    #tbfil.write('the gas below the line $\\rm DEC_{\\rm off} = \\aslope RA_{\\rm off} \\ayoff$}\n')
    tbfil.write('\\tablecomments{\\eea, \\eeb, \\eep\\ define the total $1\\sigma$ error ellipse for the FRB localization \n')
    tbfil.write('Data are taken from \\cite{Ravi19,Day20,Law20,Tendulkar17,Marcote20,Heintz2020}.} \n')
    # End
    tbfil.write('\\end{deluxetable*} \n')

    tbfil.close()
    print('Wrote {:s}'.format(outfile))


# Summary table of results
def mktab_frb_results(outfile='tab_frb_results.tex', sub=False):

    priors = [associate_defs.conservative, associate_defs.adopted]#, associate_defs.extreme]

    # Open
    tbfil = open(outfile, 'w')

    # Header
    #tbfil.write('\\clearpage\n')
    tbfil.write('\\startlongtable\n')
    tbfil.write('\\begin{deluxetable*}{cccccccccccccccc}\n')
    #tbfil.write('\\rotate\n')
    tbfil.write('\\tablewidth{0pc}\n')
    tbfil.write('\\tablecaption{Results for FRB Associations\\label{tab:results}}\n')
    tbfil.write('\\tabletypesize{\\footnotesize}\n')
    tbfil.write('\\tablehead{\\colhead{FRB}  \n')
    tbfil.write('& \\colhead{RA$_{\\rm cand}$} & \\colhead{Dec$_{\\rm cand}$} \n')
    #tbfil.write('& \\colhead{DM$_{\\rm FRB}$} 
    tbfil.write('& \\colhead{$\\theta$} \n')
    tbfil.write('& \\colhead{\\halflight}  \n')
    tbfil.write('& \\colhead{\\gmag} \n')
    tbfil.write('& \\colhead{Filter} & \\colhead{\\pchance} \n')
    tbfil.write('& \\colhead{\\PO} & \\colhead{\\POx} \n')
    tbfil.write('& \\colhead{\\PU} & \\colhead{\\PUx} \n')
    tbfil.write('\\\\')
    #tbfil.write("& (deg) & (deg) & ($''$) & & (deg) & (deg) & ($''$) & ($''$) & ($''$) & ($''$) & (mag)\n")
    #tbfil.write("\\\\ (1) & (2) & (3) & (4) & (5) & (6) & (7) & (8) & (9) & (10) & (11) & (12) & (13) & (14) & (15)")
    tbfil.write('} \n')

    tbfil.write('\\startdata \n')

    # Loop on priors
    for pp, prior in enumerate(priors):
        tbfil.write('\\cutinhead{' + '{}'.format(prior['name']) + '} \n')
        # Loop on FRBs
        for ss, frb_name in enumerate(associate_defs.frb_list):
            # Config
            config = getattr(frbs, frb_name.lower())
            config['skip_bayesian'] = True
            # Could do this outside the prior loop
            frbA = frbassociate.run_individual(config)
            # Set priors
            frbA.calc_priors(prior['U'], method=prior['O'])
            prior['theta']['ang_size'] = frbA.candidates.half_light.values
            frbA.set_theta_prior(prior['theta'])

            # Calculate p(O_i|x)
            frbA.calc_POx()
            # Reverse Sort
            frbA.candidates = frbA.candidates.sort_values('P_Ox', ascending=False)

            # Loop Galaxy
            first = True
            for cc, cand in frbA.candidates.iterrows():
                if first:
                    sline = frbA.frb.frb_name
                    first = False
                else:
                    sline = ""
                #
                sline += '& {:0.4f} & ${:0.4f}$'.format(cand.ra, cand.dec)
                sline += '& {:0.1f}'.format(cand.separation)
                sline += '& {:0.2f}'.format(cand.half_light)
                # m
                sline += '& {:0.2f}'.format(cand[frbA.filter])
                sline += '& {:s}'.format(frbA.filter.replace('_', '\_'))
                # Pchance
                sline += '& {:0.4f}'.format(cand.P_c)
                # P(M), P(M|x)
                sline += '& {:0.4f}'.format(cand.P_O)
                sline += '& {:0.4f}'.format(cand.P_Ox)
                # P(S), P(S|X)
                sline += '& {:0.4f}'.format(frbA.prior_U)
                sline += '& {:0.4f}'.format(frbA.P_Ux)

                tbfil.write(sline + '\\\\ \n')

    # End end
    tbfil.write('\\hline \n')


    tbfil.write('\\enddata \n')

    #tbfil.write('\\tablenotetext{a}{Spectroscopic redshifts are reported to 4 significant digits.  Photometric to 2.} \n')
    #tbfil.write('the gas below the line $\\rm DEC_{\\rm off} = \\aslope RA_{\\rm off} \\ayoff$}\n')
    #tbfil.write('\\tablecomments{Column 1: FRB source. Columns 2 and 3: R.A. and Decl. of the FRB (J2000). Column 4: FRB error ellipse. Column 5: FRB classication. Repeating = yes(y)/no(n). Column 6 and 7: R.A. and Dec. of the associated host galaxy (J2000). Column 8: projected angular offset of the FRB to the host galaxy center. Column 9: association radius $\delta x$ \citep{Tunnicliffe14}. Column 10: angular effective radius of the host measured from a sersic model using GALFIT \citep{galfit} on the $i$-band images (or equivalent). Column 11: effective search radius \citep{Bloom02}. Column 12: measured apparent magnitude of the host. Column 13: filter used for the magnitude measurement. Column 14: probability of chance coincidence using the \citet{Bloom02} formalism. Column 15: sample designations following the criteria outlined in $\S$~\\ref{ssec:associate}.}\n')
    # End
    tbfil.write('\\end{deluxetable*} \n')

    tbfil.close()
    print('Wrote {:s}'.format(outfile))


# Summary table of results
def mktab_sandbox(outfile='tab_sandbox.tex', sub=False):

    # Open
    tbfil = open(outfile, 'w')

    # Header
    #tbfil.write('\\clearpage\n')
    tbfil.write('\\startlongtable\n')
    tbfil.write('\\begin{deluxetable}{cccccccccccccccc}\n')
    #tbfil.write('\\rotate\n')
    tbfil.write('\\tablewidth{0pc}\n')
    tbfil.write('\\tablecaption{Sandbox Analysis\\label{tab:sandbox}}\n')
    tbfil.write('\\tabletypesize{\\footnotesize}\n')
    tbfil.write('\\tablehead{\\colhead{Sandbox}  \n')
    tbfil.write('& \\colhead{\\PO} & \\colhead{\\PU} & \\colhead{\\poffset}\n')
    tbfil.write('& \\colhead{\\thmax/\\halflight}\n')
    tbfil.write('& \\colhead{f(T+secure)} & \\colhead{TP} \n')
    tbfil.write('\\\\')
    #tbfil.write("& (deg) & (deg) & ($''$) & & (deg) & (deg) & ($''$) & ($''$) & ($''$) & ($''$) & (mag)\n")
    #tbfil.write("\\\\ (1) & (2) & (3) & (4) & (5) & (6) & (7) & (8) & (9) & (10) & (11) & (12) & (13) & (14) & (15)")
    tbfil.write('} \n')

    tbfil.write('\\startdata \n')

    aprior = associate_defs.adopted.copy()
    cprior = associate_defs.conservative.copy()

    priors = []
    SBs = []
    stat_files = []
    truth_files = []

    # Exploring P(O)
    # SB-1, inverse prior
    SBs.append('SB-1')
    prior = copy.deepcopy(aprior)
    priors.append(prior.copy())
    stat_files.append('../Analysis/SandBox/mU_oU2_1arcsec/mU_oU2_1arcsec_A_stats.csv')
    truth_files.append('../Analysis/SandBox/mU_oU2_1arcsec/frbs_with_galaxies_100000.csv')
    # SB-1, inverse1
    SBs.append('SB-1')
    prior = copy.deepcopy(aprior)
    prior['O'] = 'inverse1'
    priors.append(prior.copy())
    stat_files.append('../Analysis/SandBox/mU_oU2_1arcsec/mU_oU2_1arcsec_A1_stats.csv')
    truth_files.append('../Analysis/SandBox/mU_oU2_1arcsec/frbs_with_galaxies_100000.csv')
    # SB-1, inverse2
    SBs.append('SB-1')
    prior = copy.deepcopy(aprior)
    prior['O'] = 'inverse2'
    priors.append(prior.copy())
    stat_files.append('../Analysis/SandBox/mU_oU2_1arcsec/mU_oU2_1arcsec_A2_stats.csv')
    truth_files.append('../Analysis/SandBox/mU_oU2_1arcsec/frbs_with_galaxies_100000.csv')

    # Exploring priors
    # SB-1, conservative
    SBs.append('SB-1')
    prior = copy.deepcopy(cprior)
    priors.append(prior.copy())
    stat_files.append('../Analysis/SandBox/mU_oU2_1arcsec/mU_oU2_1arcsec_C_stats.csv')
    truth_files.append('../Analysis/SandBox/mU_oU2_1arcsec/frbs_with_galaxies_100000.csv')
    # SB-1, P(U)>0
    SBs.append('SB-1')
    prior = copy.deepcopy(aprior)
    prior['U'] = 0.05
    priors.append(prior.copy())
    stat_files.append('../Analysis/SandBox/mU_oU2_1arcsec/mU_oU2_1arcsec_A0.05_stats.csv')
    truth_files.append('../Analysis/SandBox/mU_oU2_1arcsec/frbs_with_galaxies_100000.csv')


    # SB-2
    SBs.append('SB-2')
    prior = copy.deepcopy(aprior)
    prior['theta']['method'] = 'uniform'
    prior['theta']['max'] = 2.
    priors.append(prior.copy())
    stat_files.append('../Analysis/SandBox/mag_20-23_oU2_locU_0.1-1/mag_20-23_oU2_locU_0.1-1_A_stats.csv')
    truth_files.append('../Analysis/SandBox/mag_20-23_oU2_locU_0.1-1/frbs_with_galaxies_46699_mag_20-23_oU2_locU_0.1-1.csv')
    # SB-3
    SBs.append('SB-3')
    prior = copy.deepcopy(aprior)
    prior['theta']['method'] = 'core'
    priors.append(prior.copy())
    stat_files.append('../Analysis/SandBox/mag_20-23_oCore_locU_0.1-1/mag_20-23_oCore_locU_0.1-1_A_stats.csv')
    truth_files.append('../Analysis/SandBox/mag_20-23_oCore_locU_0.1-1/frbs_with_galaxies_46699_mag_20-23_oCore_locU_0.1-1.csv')
    # SB-4
    SBs.append('SB-4')
    prior = copy.deepcopy(aprior)
    prior['theta']['method'] = 'exp'
    priors.append(prior.copy())
    stat_files.append('../Analysis/SandBox/mag_20-23_oExp_locU_0.1-1/mag_20-23_oExp_locU_0.1-1_A_stats.csv')
    truth_files.append('../Analysis/SandBox/mag_20-23_oExp_locU_0.1-1/frbs_with_galaxies_46699_mag_20-23_oExp_locU_0.1-1.csv')

    # SB-U
    SBs.append('SB-5')
    prior = copy.deepcopy(aprior)
    prior['U'] = 0.10
    priors.append(prior.copy())
    stat_files.append('../Analysis/SandBox/PU10/PU10mag25_A0.10_stats.csv')
    truth_files.append('../Analysis/SandBox/PU10/frbs_with_galaxies_50000_mag_20-25_oU2_locU_0.1-1_PU10.csv')

    # Loop on priors
    for prior, SB, stats_file, truth_file in zip(priors, SBs, stat_files, truth_files):
        if not os.path.isfile(stats_file):
            print("Skipping: {}".format(stats_file))
            continue
        sline = ''
        # SB
        sline += SB

        # P(O)
        sline += '&' + prior['O']

        # P(U)
        if prior['U'] >= 0:
            if prior['U'] == 0:
                sline += '& {}'.format(int(prior['U']))
            else:
                sline += '& {:0.2f}'.format(prior['U'])
        else:
            sline += '&' + 'chance'

        # P(theta)
        sline += '&' + prior['theta']['method']

        # theta_max
        sline += '& {}'.format(int(prior['theta']['max']))

        # Load
        truth = pandas.read_csv(truth_file, index_col=0)
        true_coords = SkyCoord(ra=truth.ra.values, dec=truth.dec.values, unit='deg')
        stats = pandas.read_csv(stats_file, index_col=0)
        pred_coords = SkyCoord(ra=stats.best_ra.values, dec=stats.best_dec.values, unit='deg')
        # Successful matches
        if SB == 'SB-5':
            iPU = stats.PU > stats.max_POx
            idx, d2d, _ = match_coordinates_sky(pred_coords, true_coords, nthneighbor=1)
            #
            correct = np.zeros(len(stats), dtype=bool)
            good_iPU = stats[iPU].iFRB > 45000
            correct[np.where(iPU)[0][good_iPU]] = True
            # Others
            good_notiPU = d2d[~iPU] < 0.01*units.arcsec
            correct[np.where(~iPU)[0][good_notiPU]] = True
        else:
            sep = pred_coords.separation(true_coords).to('arcsec')
            correct = sep < 0.01*units.arcsec
        stats['correct'] = correct

        # Max
        isecure = stats.max_POx > associate_defs.POx_secure
        if SB == 'SB-5':
            isecure = isecure | (stats.PU > associate_defs.POx_secure)
        nobj = int(np.sum(isecure))
        nhits = int(np.sum(stats.correct.values[isecure]))

        # f(True + secure)
        ncorrect = np.sum(stats.correct.values)
        Pmax = stats.max_POx[stats.correct]
        nsecure = np.sum(Pmax > associate_defs.POx_secure)
        sline += '& {:0.2f}'.format(nsecure/ncorrect)

        # TP
        true_pos = float(nhits)/nobj
        sline += '& {:0.2f}'.format(true_pos)

        # Print
        tbfil.write(sline + '\\\\ \n')

    # End end
    tbfil.write('\\hline \n')


    tbfil.write('\\enddata \n')

    #tbfil.write('\\tablenotetext{a}{Spectroscopic redshifts are reported to 4 significant digits.  Photometric to 2.} \n')
    #tbfil.write('the gas below the line $\\rm DEC_{\\rm off} = \\aslope RA_{\\rm off} \\ayoff$}\n')
    #tbfil.write('\\tablecomments{Column 1: FRB source. Columns 2 and 3: R.A. and Decl. of the FRB (J2000). Column 4: FRB error ellipse. Column 5: FRB classication. Repeating = yes(y)/no(n). Column 6 and 7: R.A. and Dec. of the associated host galaxy (J2000). Column 8: projected angular offset of the FRB to the host galaxy center. Column 9: association radius $\delta x$ \citep{Tunnicliffe14}. Column 10: angular effective radius of the host measured from a sersic model using GALFIT \citep{galfit} on the $i$-band images (or equivalent). Column 11: effective search radius \citep{Bloom02}. Column 12: measured apparent magnitude of the host. Column 13: filter used for the magnitude measurement. Column 14: probability of chance coincidence using the \citet{Bloom02} formalism. Column 15: sample designations following the criteria outlined in $\S$~\\ref{ssec:associate}.}\n')
    # End
    tbfil.write('\\end{deluxetable} \n')

    tbfil.close()
    print('Wrote {:s}'.format(outfile))



#### ########################## #########################
#### ########################## #########################
#### ########################## #########################

# Command line execution
if __name__ == '__main__':

    #mktab_priors()
    #mktab_frbs()
    #mktab_frb_results()
    mktab_sandbox()
