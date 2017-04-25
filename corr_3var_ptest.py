
# Check whether the correlation difference of two variables to a third is
# significant
# Adapted from: https://stats.stackexchange.com/questions/151219/is-one-variable-more-correlated-than-another-with-a-third-in-r

from scipy import stats
import numpy as np
import pandas as pd


def cor_diff(a, b, c):
    return np.abs(stats.pearsonr(a, b)[0] - stats.pearsonr(a, c))[0]


def standardize(arr):
    return (arr - np.nanmean(arr)) / np.nanstd(arr)


def cortest(df, a, b, c, niter=500, seed=2378942):
    '''
    Test if the difference in the correlation of a vs. b and
    a vs. c is significant using a permutation test.
    '''

    np.random.seed(seed)

    if not isinstance(df, pd.DataFrame):
        raise TypeError("df must be a pandas DataFrame.")

    # Pull out the columns and standardize
    a_vals = standardize(np.asarray(df[a]).astype(np.float))
    b_vals = standardize(np.asarray(df[b]).astype(np.float))
    c_vals = standardize(np.asarray(df[c]).astype(np.float))

    stat = cor_diff(a_vals, b_vals, c_vals)

    sim = np.empty(niter, dtype=np.float)
    for i in xrange(niter):
        # Randomly choose which indices to switch
        inds = np.random.uniform(size=a_vals.size) > 0.5

        # Switch values between b and c
        b_samps = b_vals.copy()
        b_samps[inds] = c_vals[inds]

        c_samps = c_vals.copy()
        c_samps[inds] = b_vals[inds]

        sim[i] = cor_diff(a_vals, b_samps, c_samps)

    p_value = (sim >= stat).sum() / float(niter)

    return stat, p_value, sim
