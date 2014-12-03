
'''
Script to match coords and image shapes of 2 cubes.
Specifically, this is for creating model images of single dish data for
cleaning interferometric data.
'''

import FITS_tools as ft
from astropy.io import fits
import numpy as np


def match_regrid(filename1, filename2, return_type='hdu', reappend_dim=True,
                 remove_hist=True):
    '''
    Input two fits filenames. The output will be the projection of file 1
    onto file 2
    '''

    fits1 = fits.open(filename1)
    fits2 = fits.open(filename2)

    hdr1 = fits1[0].header.copy()
    hdr2 = fits2[0].header.copy()

    if remove_hist:
        # Remove the huge CASA history
        del hdr2["HISTORY"]

    shape1 = fits1[0].data.shape
    shape2 = fits2[0].data.shape

    # We need to alter the header to make them compatible
    if len(shape1) < len(shape2):

        hdr2["NAXIS"] = len(shape1)

        del_keys = ["NAXIS", "CTYPE", "CDELT", "CRPIX", "CUNIT", "CRVAL"]

        extra_axes = \
            [posn + 1 for posn, val in enumerate(shape2[::-1]) if val == 1]

        if reappend_dim:
            deleted_keys = {}

        for ax in extra_axes:
            for del_key in del_keys:
                if reappend_dim:
                    deleted_keys[del_key+str(ax)] = hdr2[del_key+str(ax)]

                del hdr2[del_key+str(ax)]

    # Do the matching
    regrid_img = ft.hcongrid.hcongrid(fits1[0].data, fits1[0].header, hdr2)

    # Now hack the header back together!
    if reappend_dim and len(shape1) < len(shape2):
        for key in deleted_keys:
            hdr2[key] = deleted_keys[key]

        hdr2["NAXIS"] = len(shape2)

        for _ in range(len(extra_axes)):
            regrid_img = regrid_img[np.newaxis]

    # Finally, we want to take out the important portions of fits1 header
    hdr2["TELESCOP"] = hdr1["TELESCOP"]
    hdr2["DATE-OBS"] = hdr1["DATE-OBS"]
    hdr2["DATAMAX"] = hdr1["DATAMAX"]
    hdr2["DATAMIN"] = hdr1["DATAMIN"]
    hdr2["OBSERVER"] = hdr1["OBSERVER"]
    hdr2["OBJECT"] = hdr1["OBJECT"]
    hdr2["ORIGIN"] = hdr1["ORIGIN"]
    hdr2["BMAJ"] = hdr1["BMAJ"]
    hdr2["BMIN"] = hdr1["BMIN"]
    hdr2["BPA"] = hdr1["BPA"]

    return fits.PrimaryHDU(regrid_img, header=hdr2)
