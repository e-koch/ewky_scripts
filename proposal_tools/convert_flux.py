
'''
Convert from flux units
'''

import astropy.units as u
from astropy.io.fits import HDUList
import numpy as np


# Calculate area
FWHM_TO_AREA = np.pi/(4*np.log(2))


def _to_area(major, minor, unit=u.sr):
    return (major * minor * FWHM_TO_AREA).to(unit)


def convert_flux(hdu, output=u.Jy/u.beam, copy=True):
    '''
    Given a list of data and headers, convert to the specified output units.
    '''
    if copy:
        hdu = [data.copy() for data in hdu]

    # If the output unit specified has an angular part, split it off until the
    # end. The angular part comes in with the beam area.
    if len(output.to_string().split("/")) > 1:
        split_units = output.to_string().split("/")
        flux_unit = u.Unit(split_units[0])
        ang_unit = u.Unit(split_units[1])
    else:
        flux_unit = output
        ang_unit = u.sr

    if ang_unit is u.beam:
        ang_unit = u.deg**2

    for i in range(len(hdu)):
        data, header = hdu[i].data, hdu[i].header

        # Extract necessary info from the header
        bunit = header['BUNIT'].strip()
        bmaj = header['BMAJ'] * u.deg
        bmin = header['BMIN'] * u.deg
        freq = header['RESTFREQ'] * u.Hz

        area = _to_area(bmaj, bmin, unit=ang_unit)

        # Check for MJy
        MJy_flag = False
        if bunit[0] == "M":
            MJy_flag = True
            bunit = bunit[1:]

        # Figure out input flux units
        # If this doesn't work, try u.Unit('string')
        if bunit == "K":
            bunit = u.K
        elif bunit == "Jy/beam" or bunit == "Jy/bm":
            bunit = u.Jy/u.beam
        elif bunit == "Jy/sr":
            bunit = u.Jy/u.sr

        if MJy_flag:
            bunit *= 1e6

        # Convert units
        data *= bunit

        data = data.to(flux_unit,
                       equivalencies=u.brightness_temperature(area, freq)).value

        header.update('BUNIT', output.to_string())

        hdu[i].header = header
        hdu[i].data = data

    return HDUList(hdu)
