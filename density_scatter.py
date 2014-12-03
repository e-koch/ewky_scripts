
'''
Make a density/scatter plot
'''

import numpy as np
from scipy import stats
from scipy.interpolate import interp1d
import matplotlib.pyplot as p


def make_plot(x, y, gridpts=100, verbose=False):

    values = np.vstack([x, y])
    kernel = stats.gaussian_kde(values)

    xe = np.linspace(np.nanmin(x), np.nanmax(x), gridpts)
    ye = np.linspace(np.nanmin(x), np.nanmax(x), gridpts)
    xgrid,ygrid = np.meshgrid(xe,ye)
    zsurf = kernel([xgrid.ravel(),ygrid.ravel()])
    zsurf = zsurf.reshape(xgrid.shape)

    fig, ax = p.subplots()
    # Plot surface
    ax.imshow(zsurf,cmap = 'Greys',aspect='auto',extent = [xe[0],xe[-1],ye[0],ye[-1]],origin='lower')

    # Now make the contours
    zvals = np.sort(zsurf.ravel())
    cdf = np.cumsum(zvals) / np.sum(zvals)
    interpolator = interp1d(cdf,zvals)
    cvals = interpolator(1-np.array([0.9953,0.9544,0.8427]))
    ax.contour(xe,ye,zsurf,levels=cvals,colors='k')

    # Scatter plot outliers


    if verbose:
        p.show()
    else:
        return fig