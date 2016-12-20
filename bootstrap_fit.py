
from astropy.modeling import models, fitting
import numpy as np


def fit_bootstrap(model, x, y, nboot=100, plot_corner=True, noise_slices=()):
    '''
    Fit a given astropy model and use bootstrapping to estimate the
    parameter errors.

    **Bootstrapping assumes a normal distribution in the model residuals**

    '''

    fitter = fitting.LevMarLSQFitter()

    fit_model = fitter(model, x, y)

    boot_vals = np.zeros((nboot, fit_model.parameters.size))

    def errfunc(model, x, y):
        return model(x) - y

    # residuals = errfunc(model, x, y)
    if len(noise_slices) == 0:
        noisy_data = y
    else:
        noisy_data = np.hstack([y[:noise_slices[0]], y[noise_slices[1]:]])
    sigma_res = noisy_data.std()

    for i in range(nboot):

        boot_model = model.copy()

        randomDelta = np.random.normal(0., sigma_res, y.size)

        y_rand = y + randomDelta

        boot_fit = fitter(boot_model, x, y_rand)

        boot_vals[i] = boot_fit.parameters

    if plot_corner:
        import corner
        corner.corner(boot_vals,
                      quantiles=[0.16, 0.5, 0.84],
                      show_titles=False, title_args={"fontsize": 12})

    return fit_model, boot_vals
