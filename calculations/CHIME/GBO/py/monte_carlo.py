""" Monte Carlo calculations for CHIME/FRB with GBO """
import numpy as np

from scipy import interpolate

import pandas

from matplotlib import pyplot as plt
import seaborn as sns

from IPython import embed

def generate_frbs(chime_mr_tbl:str='CHIME_mr_5Jyms_150.parquet',
                  m_r_lim:float=23.3, nsample=100000,
                  debug:bool=False):

    # Load up m_r distribution
    chime_mr = pandas.read_parquet(chime_mr_tbl)
    if debug:
        chime_mr=chime_mr[0:100000]
    n_chime = len(chime_mr)

    # Generate a CDF
    srt = np.argsort(chime_mr.m_r.values)

    f_cdf = interpolate.interp1d(
        np.arange(n_chime)/(n_chime-1),
        chime_mr.m_r.values[srt])

    # Random sample
    rand = np.random.uniform(size=nsample)
    rand_mr = f_cdf(rand)

    if debug:
        plt.clf()
        ax = plt.gca()

        # Full
        cdf = np.linspace(0., 1., 10000)
        mr = f_cdf(cdf)
        ax.plot(mr, cdf, 'k-', label='Full CDF')

        # Realization
        srt2 = np.argsort(rand_mr)
        ax.step(rand_mr[srt2], 
                np.arange(nsample)/(nsample-1),
                color='r')

        ax.legend()
        plt.show()

    # Load COSMOS

    # Cut me

# Command line execution
if __name__ == '__main__':
    generate_frbs(debug=True)