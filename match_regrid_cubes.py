
'''
Script to match coords and image shapes of 2 cubes.
Specifically, this is for creating model images of single dish data for
cleaning interferometric data.
'''

import FITS_tools as ft
from astropy.io import fits
import numpy as np
from reproject import reproject
from astropy.wcs import WCS
from spectral_cube.wcs_utils import drop_axis


def match_regrid(filename1, filename2, return_type='hdu', reappend_dim=True,
                 remove_hist=True, save_output=False, save_name='new_img'):
    '''
    Input two fits filenames. The output will be the projection of file 1
    onto file 2
    '''

    fits1 = fits.open(filename1)
    fits2 = fits.open(filename2)

    hdr1 = fits1[0].header.copy()
    hdr2 = fits2[0].header.copy()

    # new_wcs = WCS(hdr2)
    # new_wcs = drop_axis(new_wcs, 3)

    # hdr2 = new_wcs.to_header()
    hdr2["CUNIT4"] = 'km/s    '
    hdr2["CRVAL4"] = -48.1391
    hdr2["CTYPE4"] = 'VELO-LSR'
    hdr2["CDELT4"] = -1.288141

    # hdr2["CDELT1"] = hdr2["CDELT1"] * (4096/256)
    # hdr2["CDELT2"] = hdr2["CDELT2"] * (4096/256)
    # hdr2["CDELT4"] = hdr2["CDELT4"] * (205/40)

    # hdr2["NAXIS1"] = 4096
    # hdr2["NAXIS2"] = 4096
    # hdr2["NAXIS4"] = 40

    # hdr2 = "SIMPLE  =                    T /Standard FITS                                   BITPIX  =                  -32 /Floating point (32 bit)                         NAXIS   =                    3                                                  NAXIS1  =                 256                                                  NAXIS2  =                 256                                                  NAXIS3  =                  205                                                  EXTEND  =                    T                                                  BSCALE  =   1.000000000000E+00                  BZERO   =   0.000000000000E+00                                                  BMAJ    =   1.755010949241E-03                                                  BMIN    =   1.572879685296E-03                                                  BPA     =   8.867265319824E+01                                                  BTYPE   = 'Intensity'                                                           OBJECT  = 'M33P1   '                                                                                                                                            BUNIT   = 'JY/BEAM '           /Brightness (pixel) unit                         EQUINOX =   2.000000000000E+03                                                  RADESYS = 'FK5     '                                                            LONPOLE =   1.800000000000E+02                                                  LATPOLE =   3.065994166667E+01                                                  PC01_01 =   1.000000000000E+00                                                  PC02_01 =   0.000000000000E+00                                                  PC03_01 =   0.000000000000E+00                                                  PC01_02 =   0.000000000000E+00                                                  PC02_02 =   1.000000000000E+00                                                  PC03_02 =   0.000000000000E+00                                                  PC01_03 =   0.000000000000E+00                                                  PC02_03 =   0.000000000000E+00                                                  PC03_03 =   1.000000000000E+00                                                  CTYPE1  = 'RA---SIN'                                                            CRVAL1  =   2.346210000000E+01                                                  CDELT1  =  -4.166666666667E-04                                                  CRPIX1  =   2.049000000000E+03                                                  CUNIT1  = 'deg     '                                                            CTYPE2  = 'DEC--SIN'                                                            CRVAL2  =   3.065994166667E+01                                                  CDELT2  =   4.166666666667E-04                                                  CRPIX2  =   2.049000000000E+03                                                  CUNIT2  = 'deg     '                                                            CTYPE3  = 'VELO-LSR'                                                            CRVAL3  =             -48.1391                                                  CDELT3  =   6.103165421963E+03                                                  CRPIX3  =   1.000000000000E+00                                                  CUNIT3  = 'km/s    '                                                            PV2_1   =   0.000000000000E+00                                                  PV2_2   =   0.000000000000E+00                                                  RESTFRQ =   1.420405751770E+09 /Rest Frequency (Hz)                             SPECSYS = 'BARYCENT'           /Spectral reference frame                        ALTRVAL =  -4.942706867471E+04 /Alternate frequency reference value             ALTRPIX =   1.000000000000E+00 /Alternate frequency reference pixel             VELREF  =                  258 /1 LSR, 2 HEL, 3 OBS, +256 Radio                 COMMENT casacore non-standard usage: 4 LSD, 5 GEO, 6 SOU, 7 GAL                 TELESCOP= 'VLA     '                                                            OBSERVER= 'unavailable'                                                         DATE-OBS= '1997-09-19T04:59:50.000002'                                          TIMESYS = 'TAI     '                                                            OBSRA   =   2.346210000000E+01                                                  OBSDEC  =   3.065994166667E+01                                                  OBSGEO-X=  -1.601185365000E+06                                                  OBSGEO-Y=  -5.041977547000E+06                                                  OBSGEO-Z=   3.554875870000E+06                                                  DATE    = '2014-12-17T20:03:37.733514' /Date FITS file was written              ORIGIN  = 'CASA 4.2.2 (prerelease r30986)'                                      END                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             "
    # hdr2 = fits.header.Header.fromstring(hdr2)

    if remove_hist:
        # Remove the huge CASA history
        del hdr2["HISTORY"]

    shape1 = fits1[0].data.shape
    shape2 = fits2[0].data.shape[:-1]

    fits2.close()

    # We need to alter the header to make them compatible
    # if len(shape1) < len(shape2):

    #     hdr2["NAXIS"] = len(shape1)

    #     del_keys = ["NAXIS", "CTYPE", "CDELT", "CRPIX", "CUNIT", "CRVAL"]

    #     extra_axes = \
    #         [posn + 1 for posn, val in enumerate(shape2[::-1]) if val == 1]

    #     if reappend_dim:
    #         deleted_keys = {}

    #     for ax in extra_axes:
    #         for del_key in del_keys:
    #             if reappend_dim:
    #                 deleted_keys[del_key+str(ax)] = hdr2[del_key+str(ax)]

    #             del hdr2[del_key+str(ax)]

    # Do the matching
    if len(shape1) == 2:
        regrid_img = ft.hcongrid.hcongrid(fits1[0].data, fits1[0].header, hdr2)
    else:
        regrid_img = ft.regrid_cube(fits1[0].data, fits1[0].header, hdr2, specaxes=(3, 3))
        # regrid_img = reproject(fits1[0], hdr2, shape_out=(205, 256, 256))[0]

    # Now hack the header back together!
    # if reappend_dim and len(shape1) < len(shape2):
    #     for key in deleted_keys:
    #         hdr2[key] = deleted_keys[key]

    #     hdr2["NAXIS"] = len(shape2)

    #     for _ in range(len(extra_axes)):
    #         regrid_img = regrid_img[np.newaxis]

    # Finally, we want to take out the important portions of fits1 header
    # hdr2["TELESCOPE"] = hdr1["TELESCOPE"]
    # hdr2["DATE-OBS"] = hdr1["DATE-OBS"]
    # hdr2["DATAMAX"] = hdr1["DATAMAX"]
    # hdr2["DATAMIN"] = hdr1["DATAMIN"]
    # hdr2["OBSERVER"] = hdr1["OBSERVER"]
    # hdr2["OBJECT"] = hdr1["OBJECT"]
    # hdr2["ORIGIN"] = hdr1["ORIGIN"]
    # hdr2["BMAJ"] = hdr1["BMAJ"]
    # hdr2["BMIN"] = hdr1["BMIN"]
    # hdr2["BPA"] = hdr1["BPA"]

    if save_output:
        hdu = fits.PrimaryHDU(regrid_img, header=hdr2)
        hdu.writeto(save_name+".fits")

    else:
        return fits.PrimaryHDU(regrid_img, header=hdr2)


if __name__ == '__main__':

    import sys

    file1 = str(sys.argv[1])
    file2 = str(sys.argv[2])

    match_regrid(file1, file2, save_output=True, save_name='M33_vla_mask.fits')
