
from astropy.io import fits
import matplotlib.pyplot as p
from astropy.visualization import SqrtStretch, AsinhStretch
from astropy.visualization.mpl_normalize import ImageNormalize
import numpy as np
from reproject import reproject_interp
import astropy.wcs as wcs

'''
NGC 5055 over multiple tracers/wavelengths
'''

# Load in the data into a dict


Nrows = 2
Ncols = 3

fig, ax = p.subplots(Nrows, Ncols,
                     sharex=True,
                     sharey=True, num=1)

data = {}

data["HI"] = fits.open("/media/eric/Data_3/VLA/THINGS/NGC_5055/NGC_5055_NA_MOM0_THINGS.FITS")[0]
data["CO(2-1)"] = fits.open("/media/eric/Data_3/HERACLES/ngc5055_heracles_mom0.fits")[0]
data[r"H$\alpha$"] = fits.open("/media/eric/Data_3/NGC_5055/ngc5055_HA_SUB_dr4.fits")[0]
data["24 $\mu$m"] = fits.open("/media/eric/Data_3/NGC_5055/ngc5055_mips24_crop_v5-0.fits")[0]
data["250 $\mu$m"] = fits.open("/media/eric/Data_3/NGC_5055/NGC5055_kingfish_spire250_v3-0_scan.fits")[0]
data["UV"] = fits.open("/media/eric/Data_3/NGC_5055/NGA_NGC5055-fd-int.fits")[0]

# Redundant axes
hi_wcs_corr = wcs.WCS(data["HI"].header).dropaxis(3).dropaxis(2)

# Slice the HI so the centre is more visible
slices = (slice(330, 680), slice(250, 820))

ordering = ["HI", "CO(2-1)", r"H$\alpha$", "24 $\mu$m", "250 $\mu$m", "UV"]

for ii, key in enumerate(ordering):

    if key != "HI":
        # Why so many redundant axes!!
        if key == "CO(2-1)":
            arr = data[key].data.squeeze()
            mywcs = wcs.WCS(data[key].header).dropaxis(2)
            input_data = (arr, mywcs)
        else:
            input_data = data[key]
        # Regrid to the HERACLES data
        arr = reproject_interp(input_data, hi_wcs_corr,
                               shape_out=data["HI"].data.squeeze().shape)[0]
    else:
        arr = data[key].data.squeeze()

    y, x = np.unravel_index(ii, (Nrows, Ncols))
    # if key == r"H$\alpha$":
    #     vmin = 0
    #     vmax = 10
    # else:
    vmin = None
    vmax = None
    im = ax[y, x].imshow(np.arctan(arr[slices] / np.nanpercentile(arr[slices], 85)),
                         origin='lower', cmap=p.cm.gray_r,
                         vmin=vmin, vmax=vmax)

    ax[y, x].annotate(key, (0.7, 0.9),
                      xycoords='axes fraction', color='k',
                      fontsize=15.5)

# Only keep certain ticks/labels
# See https://github.com/keflavich/paper_w51_evla/blob/master/plot_codes/h77a_layers.py
for i in range(Nrows):
    for j in range(Ncols):
        if i == 0:
            ax[i, j].xaxis.set_ticks_position('top')
            p.setp(ax[i, j].get_xticklabels(), visible=False)
            ax[i, j].xaxis.set_ticklabels([])
        elif i == Nrows - 1:
            ax[i, j].xaxis.set_ticks_position('bottom')
            p.setp(ax[i, j].get_xticklabels(), visible=True)
        else:
            ax[i, j].xaxis.set_ticks_position('none')
            p.setp(ax[i, j].get_xticklabels(), visible=False)
            ax[i, j].xaxis.set_ticklabels([])

        if j == 0:
            ax[i, j].yaxis.set_ticks_position('left')
        elif j == Ncols - 1:
            ax[i, j].yaxis.set_ticks_position('right')
            p.setp(ax[i, j].get_yticklabels(), visible=False)
            ax[i, j].yaxis.set_ticklabels([])
        else:
            ax[i, j].yaxis.set_ticks_position('none')
            p.setp(ax[i, j].get_yticklabels(), visible=False)
            ax[i, j].yaxis.set_ticklabels([])

p.subplots_adjust(hspace=0,
                  wspace=0)

raw_input("Next plot?")
p.clf()

# Hi moment 0 and moment1

fig, ax = p.subplots(1, 2,
                     sharex=True,
                     sharey=True, num=1)

ax[0].imshow(data["HI"].data.squeeze(), origin='lower', cmap='afmhot')
ax[0].yaxis.set_ticks_position('none')
p.setp(ax[0].get_yticklabels(), visible=False)
ax[0].yaxis.set_ticklabels([])
ax[0].xaxis.set_ticks_position('none')
p.setp(ax[0].get_xticklabels(), visible=False)
ax[0].xaxis.set_ticklabels([])
# ax[0].annotate("Moment 0", (0.7, 0.9),
#                xycoords='axes fraction', color='k',
#                fontsize=15.5)

mom1 = fits.open("/media/eric/Data_3/VLA/THINGS/NGC_5055/NGC_5055_NA_MOM1_THINGS.FITS")[0]
ax[1].imshow(mom1.data.squeeze(), origin='lower', cmap='seismic')
ax[1].yaxis.set_ticks_position('none')
p.setp(ax[1].get_yticklabels(), visible=False)
ax[1].yaxis.set_ticklabels([])
ax[1].xaxis.set_ticks_position('none')
p.setp(ax[1].get_xticklabels(), visible=False)
ax[1].xaxis.set_ticklabels([])
# ax[1].annotate("Moment 1", (0.7, 0.9),
#                xycoords='axes fraction', color='k',
#                fontsize=15.5)

p.tight_layout()
