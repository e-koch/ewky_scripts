
'''
Find the amount of emission recovered in a high-resolution map, based on a low-resolution map
'''

from spectral_cube import SpectralCube
import numpy as np
from radio_beam import Beam
from astropy.convolution import convolve
from astropy import units as u
from astropy import wcs


class MultiResObs(object):
    """
    Object to hold observations of the same object at different resolutions.
    This is intended for matching and analyzing high-res (interferometric)
    and low-res (single dish) observations.
    """
    def __init__(self, highres, lowres):
        super(MultiResObs, self).__init__()
        self.highres = SpectralCube.read(highres)
        self.lowres = SpectralCube.read(lowres)

        self.lowbeam = self.lowres.beam
        self.highbeam = self.highres.beam

        self.combined_beam = self.lowbeam.convolve(self.highbeam)

    def match_coords(self):
        '''
        Match the spatial and spectral coordinates of the cubes.
        '''

        # Are either of the cubes in the correct frame?
        if self.highres.header['CTYPE1'] != self.lowres.header['CTYPE1']:
            raise TypeError("ctypes do not match. Are observations in the "
                            "same projection? Highres: " +
                            self.highres.header['CTYPE1'] +
                            " Lowres: " + self.lowres.header['CTYPE1'])
        elif self.highres.header['CTYPE2'] != self.lowres.header['CTYPE2']:
            raise TypeError("ctypes do not match. Are observations in the "
                            "same projection? Highres: " +
                            self.highres.header['CTYPE2'] +
                            " Lowres: " + self.lowres.header['CTYPE2'])

        # Determine which cube should be slice down so both have the same
        # spatial coverage.

        all_extrema = []

        long_slice = self.lowres.longitude_extrema > self.highres.longitude_extrema

        if long_slice[0]:
            all_extrema.append(self.highres.longitude_extrema[0])
        else:
            all_extrema.append(self.lowres.longitude_extrema[0])

        if long_slice[1]:
            all_extrema.append(self.highres.longitude_extrema[1])
        else:
            all_extrema.append(self.lowres.longitude_extrema[1])

        lat_slice = self.lowres.latitude_extrema > self.highres.latitude_extrema

        if lat_slice[0]:
            all_extrema.append(self.highres.latitude_extrema[1])
        else:
            all_extrema.append(self.lowres.latitude_extrema[1])

        if lat_slice[1]:
            all_extrema.append(self.highres.latitude_extrema[0])
        else:
            all_extrema.append(self.lowres.latitude_extrema[0])

        # Apply common slicing to both

        self.highres = \
            self.highres.subcube(xlo=all_extrema[0], xhi=all_extrema[1],
                                 ylo=all_extrema[2], yhi=all_extrema[3])

        self.lowres = \
            self.lowres.subcube(xlo=all_extrema[0], xhi=all_extrema[1],
                                ylo=all_extrema[2], yhi=all_extrema[3])

        # Now match the spectral extends

        low_spec = \
            self.highres.spectral_extrema[0] if \
            self.highres.spectral_extrema[0] > self.lowres.spectral_extrema[0]\
            else self.lowres.spectral_extrema[0]

        high_spec = \
            self.highres.spectral_extrema[1] if \
            self.highres.spectral_extrema[1] < self.lowres.spectral_extrema[1]\
            else self.lowres.spectral_extrema[1]

        self.highres = self.highres.spectral_slab(low_spec, high_spec)
        self.lowres = self.lowres.spectral_slab(low_spec, high_spec)

        return self

    def convert_to(self, unit=u.K, freq=1420.40575177*u.MHz):
        '''
        Convert both sets to common brightness units.
        '''

        convert_high = unit == self.highres.unit
        convert_low = unit == self.lowres.unit

        high_unit = self.highres.unit
        low_unit = self.lowres.unit

        if unit == u.K:

            if convert_high and 'Jy' in high_unit.name:

                self.highres = \
                    self.highres.to(unit, self.highbeam.jtok_equiv(freq))

            if convert_low and 'Jy' in low_unit.name:

                self.lowres = \
                    self.lowres.to(unit, self.lowbeam.jtok_equiv(freq))
        else:
            raise NotImplementedError("Only supporting Jy/beam -> K right now.")

    def apply_mask(self, highres_mask=None, lowres_mask=None):
        '''
        Apply a pre-made mask to either of the data cubes.
        '''

        if highres_mask is not None:
            self.highres = self.highres.with_mask(highres_mask)

        if lowres_mask is not None:
            self.lowres = self.lowres.with_mask(lowres_mask)

    def convolve_to_common(self):
        '''
        Convolve cubes to a common resolution using the combined beam.
        '''
        pass


    def flux_recovered(self):
        pass
