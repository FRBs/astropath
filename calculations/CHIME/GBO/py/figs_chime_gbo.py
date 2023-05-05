""" Figures related to CHIME/GBO calculations """
# imports
import numpy as np

import pandas

import seaborn as sns
from matplotlib import pyplot as plt


def simple_stats(): 
    pass

# Command line execution
if __name__ == '__main__':

    # Generate FRBs
    #generate_frbs('frb_monte_carlo.csv',
    #    debug=True, plots=False, nsample=10000)

    # Monte Carlo
    #run_mc('first_try.csv', debug=True)
    run_mc('frb_monte_carlo_15x3.csv', 'PATH_15x3.csv')#, debug=True)