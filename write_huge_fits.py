
'''
Save massive FITS files.
'''

from astropy.io import fits
from astropy.utils.console import ProgressBar
import numpy as np
import os


def create_huge_fits(filename, header, shape=None):

    header.tofile(filename)

    if shape is None:
        try:
            shape = tuple(header['NAXIS{0}'.format(ii)] for ii in
                          range(1, header['NAXIS'] + 1))
        except KeyError:
            raise KeyError("header does not contain the NAXIS keywords. Need "
                           "to provide a shape.")
    with open(filename, 'rb+') as fobj:
        fobj.seek(len(header.tostring()) + (np.product(shape) *
                                            np.abs(header['BITPIX'] // 8)) - 1)
        fobj.write(b'\0')


def write_huge_fits(cube, filename, verbose=True, dtype=">f4"):
    '''
    Write a SpectralCube to a massive FITS file. Currently only 3D objects
    will work. A more general approach where the longest axis is iterated over
    is probably a good idea.

    cube : np.ndarray or SpectralCube
        Data to be written.

    filename : str
        Name of the FITS file to write to. If it does not already exist,
        `create_huge_fits` is called first to create an empty FITS file.
    '''

    if not os.path.exists(filename):
        create_huge_fits(filename, cube.header)

    # Open the FITS and write the values out channel-by-channel
    hdu = fits.open(filename, mode='update')

    nchans = cube.shape[0]

    if verbose:
        iterat = ProgressBar(nchans)
    else:
        iterat = xrange(nchans)

    for i in iterat:
        hdu[0].data[i, :, :] = cube[i, :, :].value.astype(dtype)
        hdu.flush()

    # And write the header
    hdu[0].header = cube.header
    hdu.flush()

    hdu.close()
