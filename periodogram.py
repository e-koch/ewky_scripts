
'''
Script mostly written by @astrolebs.
'''

# IMPORT STUFFS
import numpy as np
import matplotlib.pyplot as p
import numpy.random as ra

from scipy import linalg


#########################################################
# CREATE BAYESPERIOD FUNCTION

def BayesPeriod(inputdata, time, periods):

    # DEFINE CONSTANTS
    N = len(time)
    data = inputdata - np.mean(inputdata)
    T = time

    # FIT DATA DEPENDENT ON TIME IN A POLYNOMIAL OF ONE DEGREE
    mod = np.polyfit(T, data, 1)

    # CALCULATE EXPECTED VALUES
    datanew = data - (mod[0] * T + mod[1])

    sigma = np.var(data)
    DataPower = sum(data ** 2) / N
    freqs = 1 / periods
    lnprob = []
    z = 1

    # CALCULATE lnprob VECTOR
    for num in range(0, len(freqs) - 1):
        theta = (0.5) * np.arctan2(sum(np.sin(4 * np.pi *
                                              freqs[num] * time)) * z, sum(np.cos(4 * np.pi * freqs[num] * time) * z))
        costerm = sum(data * np.cos(2 * np.pi * freqs[num] * time - theta))
        sinterm = sum(data * np.sin(2 * np.pi * freqs[num] * time - theta))
        Cterm = sum(np.cos(2 * np.pi * freqs[num] * time - theta) ** 2)
        Sterm = sum(np.sin(2 * np.pi * freqs[num] * time - theta) ** 2)
        Hsquare = (costerm ** 2) / Cterm + (sinterm ** 2) / Sterm
        lnprob[num] = (-0.5) * np.log(Cterm) - (0.5) * \
            np.log(Sterm) + np.log(N * DataPower - Hsquare) * (1 - N / 2)

    lnprob = lnprob / sigma

    return lnprob

#######################################################
# CREATE PERIOD FUNCTION


def Period(data, time, returnMax=True, plotPeriod=False, realPeriod=0):

    # GET TestPeriod AND lnprob
    TestPeriod = np.exp(np.linspace(0, 10, num=10000, endpoint=True))
    lnprob = BaysPeriod(data, time, TestPeriod)

    # CALCULATE MAX OF lnprob FOR CONVENIENCE AND MY EASE OF MIND
    lpmax = max(lnprob)

    # PLOTS
    if plotPeriod:
        # PLOT PRODUCTES plot(x = log10(TestPeriod), y = (lnprob-max(lnprob)))
        p.plot(np.log10(TestPeriod), (lnprob - lpmax))
        p.title("Bayes Period")
        if np.not_equal(realPeriod, 0):
            # PLOT VERTICAL LINE AT VALUE = realPeriod
            '''I AM NOT SURE IF THIS IS THE RIGHT COMMAND '''
            p.vlines(realPeriod)
        if returnMax:
            # test period[which.max(lnprob)]
            '''AGAIN, NOT SURE THAT THIS IS THE RITHGT COMMAND :) '''
            print(TestPeriod[lnprob.index(lpmax)])
        p.show()

# if __name__=="__main__":
# 	import sys
# 	fib(int(sys.argv[1]))
