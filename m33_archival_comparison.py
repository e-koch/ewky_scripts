
'''
Create comparison figure of archival VLA M33 data
'''

from spectral_cube import SpectralCube
from signal_id import Noise, RadioMask
import aplpy
from matplotlib import pyplot as p
from astropy import units as u

# Load in the cubes
cube = SpectralCube.read("M33_206_b_c_HI.fits")
cube = cube.with_mask(cube != 0*u.Jy)

old_cube = SpectralCube.read("../old_imaging/m33_hi.fits")

# Make a comparison of single channels
# Using channel 100 in new cube

chan_100 = cube.spectral_slab(cube.spectral_axis[100], cube.spectral_axis[100])

old_chan_100 = \
    old_cube.spectral_slab(cube.spectral_axis[100], cube.spectral_axis[100])


chan_fig = p.figure()

sub1 = aplpy.FITSFigure(chan_100.hdu, subplot=(1, 2, 1), figure=chan_fig)
sub2 = aplpy.FITSFigure(old_chan_100.hdu, subplot=(1, 2, 2), figure=chan_fig)

sub1.show_grayscale(invert=True)
sub2.show_grayscale(invert=True)

sub1.show_colorbar()
sub2.show_colorbar()

sub1.colorbar.set_axis_label_text("Intensity (Jy/beam)")
sub2.hide_yaxis_label()

chan_fig.tight_layout()

raw_input("Continue?")