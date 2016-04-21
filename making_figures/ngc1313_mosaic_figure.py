
from aplpy import FITSFigure, make_rgb_image
import os

'''
For proposal: 2016.1.00750.S
'''

path = "/media/eric/Data_3/NGC1313/"
fits_cube = os.path.join(path, "ngc1313_noao.fits")
fits_twod = os.path.join(path, "ngc1313_noao_2d.fits")
rgb_image = os.path.join(path, "ngc1313_noao_rgb.png")

# Make RGB image if it isn't there
if not os.path.exists(rgb_image):
    make_rgb_image(fits_cube, rgb_image)

fig = FITSFigure(fits_twod)
fig.show_rgb(rgb_image)
# fig.show_grayscale(invert=True)
# ALMA regions
fig.show_regions("ngc1313_targets.reg")
# Hubble ACS field outlines
fig.show_regions(os.path.join(path, "HST/final_mosaic/coverage_outline.reg"))
fig.hide_tick_labels()
fig.hide_xaxis_label()
fig.hide_yaxis_label()

# Should be zoomed in on before saving
# fig.save(os.path.join(path, "ngc1313_noao_field_overlays.png"))

# Also needs to be re-orientated such that N is actually up-sh
