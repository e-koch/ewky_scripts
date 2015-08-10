
'''
Create comparison figure of archival VLA M33 data
'''

from spectral_cube import SpectralCube
from signal_id import Noise, RadioMask
from radio_beam import Beam
import aplpy
from matplotlib import pyplot as p
from astropy import units as u
from astropy.io import fits

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

# Now create good masks and derive 0th moments.


new_beam = Beam.from_fits_header(cube.header)
noise_cube = Noise(cube, beam=new_beam)

new_noise = noise_cube.scale

cube_masked = cube.with_mask(cube > new_noise*u.Jy)

old_beam = Beam.from_fits_header(old_cube.header)
old_noise_cube = Noise(old_cube, beam=old_beam)

old_noise = old_noise_cube.scale

old_cube_masked = old_cube.with_mask(old_cube > old_noise*u.Jy)

# Load in the broad clean mask used
clean_mask = fits.getdata("../../../Arecibo/M33_newmask.fits")

# Need to match the dims
clean_mask = clean_mask.squeeze()
clean_mask = clean_mask[11:195, ::-1, ::-1]
clean_mask = clean_mask[:, 595:3504, 1065:3033]

mask = RadioMask(cube_masked)
mask.intersection(clean_mask)
# mask.remove_small_regions()
mask.open(iterations=3)
mask.close(iterations=3)
mask.dilate(iterations=6)

 old_arecibo_mask = fits.getdata("../old_imaging/M33_arecibo_mask_old_AT0206.fits")

old_mask = RadioMask(old_cube_masked)
old_mask.intersection(old_arecibo_mask)
old_mask.remove_small_regions()
