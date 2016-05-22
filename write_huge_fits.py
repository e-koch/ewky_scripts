
'''
Save massive FITS files.
'''

from astropy.io import fits
from astropy.utils.console import ProgressBar
import numpy as np
import functools
import operator
import os


def create_huge_fits(shape, filename, header=None, nblocks=4,
                     dtype=float):
    '''
    Creates a massive empty FITS file that can be then written to
    in slice (or something that doesn't require reading it all in).

    Code from the handy astropy FAQ:
    http://astropy.readthedocs.org/en/latest/io/fits/appendix/faq.html#how-can-i-create-a-very-large-fits-file-from-scratch

    Parameters
    ----------
    shape : tuple
        Shape of the data.
    filename : str
        Outputted FITS file name.
    header: FITS header, optional
        Provide a header for the FITS file. If given, the shape in the header
        will be used (if it is given) before the given shape is used.
    '''

    if header is not None:
        try:
            naxis = header["NAXIS"]

            shape = []

            for axis in range(naxis, 0, -1):
                shape.append(header["NAXIS" + str(axis)])

            shape = tuple(shape)
        except KeyError:
            Warning("Cannot extract shape info from header. Using given"
                    " shape instead.")

    naxis = len(shape)

    # Create an array with the right number of dimensions.

    inp_data = np.zeros((100, ) * naxis, dtype=dtype)

    # Make hdu and pad header with enough header blocks
    hdu = fits.PrimaryHDU(data=inp_data)

    hdr = hdu.header

    while len(hdr) < (36 * nblocks - 1):
        hdr.append()

    # Set the actual shape in the header
    hdr["NAXIS"] = naxis

    for axis in range(1, naxis + 1):
        hdr["NAXIS" + str(axis)] = shape[naxis - axis]

    # Save the header
    filename = os.path.splitext(filename)[0] + ".fits"
    hdr.tofile(filename)

    # Now stream some zeros into it!
    nelements = functools.reduce(operator.mul, shape, 1)
    with open(filename, 'rb+') as fobj:
        fobj.seek(len(hdr.tostring()) + (8 * nelements) - 1)
        fobj.write('\0')


def write_huge_fits(cube, filename, verbose=True):
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
        create_huge_fits(cube.shape, filename, header=cube.header)

    # Open the FITS and write the values out channel-by-channel
    hdu = fits.open(filename, mode='update')

    nchans = cube.shape[0]

    if verbose:
        iterat = ProgressBar(nchans)
    else:
        iterat = xrange(nchans)

    for i in iterat:
        hdu[0].data[i, :, :] = cube[i, :, :]
        hdu.flush()

    # And write the header
    hdu[0].header = cube.header
    hdu.flush()

    hdu.close()
