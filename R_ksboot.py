
'''
Function to run KS bootstrap test from R
'''

from rpy2.robjects.packages import importr
from rpy2.robjects.numpy2ri import numpy2ri as np2r

# Load it in
matching = importr("Matching")
base = importr("base")

def ks_boot(sample1, sample2, nboots=20000, alternative="two.sided", verbose=False):

    # This apparently converts numpy arrays to R objects

    print_steps = 0
    if verbose:
        print_steps = 2

    sample1 = np2r(sample1)
    sample2 = np2r(sample2)

    ks = matching.ks_boot(sample1, sample2, nboots, alternative, print_steps)

    if verbose:
        print(base.summary(ks))

    pval = float(ks[0][0])
    statistic = float(ks[1].rx2("statistic")[0])

    return statistic, pval

