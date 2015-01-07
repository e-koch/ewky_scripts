
'''
Convert Jy/beam to Jy/sr
'''

from astropy import units as u
from astropy.io import fits
import numpy as np


def convert_angular(filename, beam):

    hdu = fits.open(filename)

    # Assume -1 extension is the image and 2 is the error

    signal_map = hdu[-1].copy()
    error_map = hdu[2].copy()

    primary = hdu[0].copy()

    # Only need this for Jy/bm right now
    beam = beam * u.arcsec

    # Convert beam to sr and divide by 1e6 to get MJy/sr
    signal_map.data = convert_func(signal_map.data, beam).value / 1e6
    error_map.data = convert_func(error_map.data, beam).value / 1e6

    # Update the headers
    signal_map.header.update("BUNIT", "MJy/sr")
    error_map.header.update("BUNIT", "MJy/sr")

    primary.header.update("BUNIT", "MJy/sr")
    del primary.header["COMMENT"]
    primary.header["COMMENT"] = "The file contains 2 extensions:"
    primary.header["COMMENT"] = "1 - signal map (weighted, extension 4 in original)"
    primary.header["COMMENT"] = "2 - error map (extension 1 in original)"

    # Remove .fits
    filename = filename[:-5]

    # Save the signal and error maps
    new_hdu = fits.HDUList([primary, signal_map, error_map])

    new_hdu.writeto(filename + "_converted.fits")


def convert_func(arr, beam):
    factor = (2*np.pi / (8*np.log(2))) * (beam**2).to(u.sr)

    return arr / factor


if __name__ == "__main__":

    import sys

    filename = sys.argv[1]

    # Assuming arcsec here
    beam = float(sys.argv[2])

    convert_angular(filename, beam)
