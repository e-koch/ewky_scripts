
from astropy.io import fits
import sys

'''
Given a FITS cube, spit out each channel as its own file. Doesn't handle
proper changing of the spectral dimension in the header, because I don't
need it ATM.
'''

cube_name = sys.argv[1]
output_dir = sys.argv[2]
spec_axis = int(sys.argv[3])

try:
    hdu_posn = int(sys.argv[4])
except IndexError:
    hdu_posn = 0

if output_dir[-1] != "/":
    output_dir += "/"

cube = fits.open(cube_name)[hdu_posn]

chan_slice = [slice(None)] * len(cube.data.shape)

hdr_spec_axis = len(cube.data.shape) - spec_axis - 1

for chan in range(1, cube.shape[spec_axis]+1):

    print("Saving channel %s of %s" % (chan, cube.shape[spec_axis]))

    chan_hdr = cube.header.copy()

    chan_hdr['NAXIS'+str(hdr_spec_axis)] = 1

    chan_slice[spec_axis] = slice(chan-1, chan, None)

    chan_hdu = fits.PrimaryHDU(cube.data[chan_slice], header=chan_hdr)

    chan_hdu.writeto(output_dir+cube_name[:-5]+"_channel_"+str(chan)+".fits")
