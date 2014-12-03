
'''
Fit a broken power law with a background (noise)
'''

import numpy as np
import matplotlib.pyplot as p


def bkplw(k, A, brk, alp1, alp2):
    '''
    Adapted from https://github.com/keflavich/plfit.
    '''
    A = 10**A

    scale2loc = np.argmin(np.abs(k - brk))
    A2 = A * k[scale2loc]**(alp1-alp2)
    return np.log10(A * k**alp1 * (k > brk) + A2 * k**alp2 * (k <= brk))


def fit_func(pars, data, k):
    return -np.nansum((data - bkplw(k, *pars))**2.)

if __name__ == "__main__":

    # Load in a VCS set

    import astropy.io.fits as fits

    cube, hdr = \
        fits.getdata("/Users/ekoch/Dropbox/AstroStatistics/ngc1333.13co.fits",
                     header=True)

    from turbustat.statistics import VCS

    vcs = VCS(cube, hdr, phys_units=False)

    vcs.run(verbose=False)

    k = vcs.vel_freqs
    plaw = np.log10(vcs.ps1D)
    plaw = plaw[np.nonzero(k)]
    k = k[np.nonzero(k)]

    from scipy.optimize import curve_fit

    results = curve_fit(bkplw, k, plaw, p0=(np.log10(1e3), 0.1, 0.0, -2.),
                        maxfev=100000)
    print results[0]
    p.loglog(vcs.vel_freqs, vcs.ps1D, 'bD')
    p.loglog(k, 10**bkplw(k, *results[0]), 'r')
    p.show()

    from emcee import EnsembleSampler

    ndim = 4
    nwalkers = 100

    pos = [results[0] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

    sampler = EnsembleSampler(nwalkers, ndim, fit_func, args=(plaw, k))
    sampler.run_mcmc(pos, 1e3)

    samples = sampler.chain.reshape((-1, ndim))
    # Remove the burn-in
    samples = samples[200:, :]

    import triangle

    fig = triangle.corner(samples)
    p.show()

    print samples.mean(axis=0)
    p.loglog(vcs.vel_freqs, vcs.ps1D, 'bD')
    p.loglog(k, 10**bkplw(k, *results[0]), 'r')
    p.loglog(k, 10**bkplw(k, *samples.mean(axis=0)), 'g')
    p.grid(True)
    p.show()
