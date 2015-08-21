
'''
Find the amount of emission recovered in a high-resolution map, based on a low-resolution map
'''

from spectral_cube import SpectralCube
import numpy as np
from radio_beam import Beam
from astropy.convolution import convolve_fft
from astropy import units as u
from astropy import wcs
from FITS_tools.header_tools import wcs_to_platescale
import matplotlib.pyplot as p
from scipy.interpolate import interp1d
import dask.array as da


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

        self.highres_convolved = None
        self.lowres_convolved = None

        self.lowbeam = self.lowres.beam
        self.highbeam = self.highres.beam

        self.combined_beam = self.lowbeam.convolve(self.highbeam)

    def apply_mask(self, highres_mask=None, lowres_mask=None):
        '''
        Apply a pre-made mask to either of the data cubes.
        '''

        if highres_mask is not None:
            self.highres = self.highres.with_mask(highres_mask)

        if lowres_mask is not None:
            self.lowres = self.lowres.with_mask(lowres_mask)

    def match_coords(self):
        '''
        Match the spatial and spectral coordinates of the cubes.
        '''

        # TODO: use Skycoords to convert extrema into the same frame

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

        limits_high = []
        limits_low = []

        low_long = self.lowres.longitude_extrema
        high_long = self.highres.longitude_extrema

        low_lat = self.lowres.latitude_extrema
        high_lat = self.highres.latitude_extrema

        if low_long[0] < high_long[0]:
            limits_low.append(high_long[0])
            limits_high.append("min")
        else:
            limits_high.append(low_long[0])
            limits_low.append("min")

        if low_long[1] > high_long[1]:
            limits_low.append(high_long[1])
            limits_high.append("max")
        else:
            limits_high.append(low_long[0])
            limits_low.append("min")

        if low_lat[1] > high_lat[1]:
            limits_low.append(high_lat[1])
            limits_high.append("min")
        else:
            limits_high.append(low_lat[0])
            limits_low.append("min")

        if low_lat[0] < high_lat[0]:
            limits_low.append(high_lat[0])
            limits_high.append("max")
        else:
            limits_high.append(low_lat[0])
            limits_low.append("max")

        # Apply common slicing to both

        self.highres = \
            self.highres.subcube(xlo=limits_high[0], xhi=limits_high[1],
                                 ylo=limits_high[2], yhi=limits_high[3])

        self.lowres = \
            self.lowres.subcube(xlo=limits_low[0], xhi=limits_low[1],
                                ylo=limits_low[2], yhi=limits_low[3])

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

    def convert_to(self, unit=u.K, freq=1420.40575177*u.MHz):
        '''
        Convert both sets to common brightness units.
        '''

        convert_high = unit != self.highres.unit
        convert_low = unit != self.lowres.unit

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

    def convolve_to_common(self, verbose=False, use_dask=True,
                           block=(256, 256)):
        '''
        Convolve cubes to a common resolution using the combined beam.
        '''

        # Create convolution kernels from the combined beam
        conv_kernel_high = \
            self.combined_beam.as_kernel(wcs_to_platescale(self.highres.wcs))

        # if use_dask:
        #     assert np.alltrue([bl > kern for bl, kern in
        #                        zip(block, conv_kernel_high.shape)])

        conv_kernel_low = \
            self.combined_beam.as_kernel(wcs_to_platescale(self.lowres.wcs))

        # if use_dask:
        #     assert np.alltrue([bl > kern for bl, kern in
        #                        zip(block, conv_kernel_low.shape)])

        high_pad = np.ceil(conv_kernel_high.shape[0] / 2).astype(int)

        highres_convolved = np.empty(self.highres.shape)

        high_chans = len(self.highres.spectral_axis)

        if verbose:
            print("Convolving high resolution cube.")

        if use_dask:
            highres_convolved = auto_dask_map(self.highres, blocks=block,
                                              args=[conv_kernel_high])

        # for chan in range(high_chans):
        #     if verbose:
        #         print("On Channel: "+str(chan)+" of "+str(high_chans))

        #     if use_dask:
        #         da_arr = \
        #             da.from_array(np.pad(self.highres.filled_data[chan, :, :],
        #                                  high_pad, padwithnans),
        #                           chunks=block)

        #         highres_convolved[chan, high_pad:-high_pad,
        #                           high_pad:-high_pad] = \
        #             da_arr.map_overlap(
        #                 lambda a:
        #                     convolve_fft(a,
        #                                  conv_kernel_high,
        #                                  boundary='fill',
        #                                  interpolate_nan=True,
        #                                  normalize_kernel=True),
        #                 depth=2*high_pad,
        #                 boundary=np.nan).compute()

        #     else:
        #         highres_convolved[chan, :, :] = \
        #             convolve_fft(self.highres.filled_data[chan, :, :],
        #                          conv_kernel_high, boundary='fill',
        #                          interpolate_nan=True, normalize_kernel=True)
        update_high_hdr = \
            _update_beam_in_hdr(self.highres.header, self.combined_beam)

        self.highres_convolved = \
            SpectralCube(highres_convolved*self.highres.unit, self.highres.wcs,
                         header=update_high_hdr)

        # Cleanup a bit

        del highres_convolved

        # Now the low resolution data

        lowres_convolved = np.empty(self.lowres.shape)

        low_chans = len(self.lowres.spectral_axis)

        if verbose:
            print("Convolving low resolution cube.")

        if use_dask:
            lowres_convolved = auto_dask_map(self.highres, blocks=block,
                                             args=[conv_kernel_low])

        # for chan in range(low_chans):
        #     if verbose:
        #         print("On Channel: "+str(chan)+" of "+str(low_chans))

        #     lowres_convolved[chan, :, :] = \
        #         convolve_fft(self.lowres.filled_data[chan, :, :],
        #                      conv_kernel_low, boundary='fill',
        #                      interpolate_nan=True, normalize_kernel=True)

        update_low_hdr = \
            _update_beam_in_hdr(self.lowres.header, self.combined_beam)

        self.lowres_convolved = \
            SpectralCube(lowres_convolved*self.lowres.unit, self.lowres.wcs,
                         header=update_low_hdr)

    def flux_recovered(self, plot=True, filename=None, enable_interp=True,
                       interp_to='lower'):
        '''
        Check the amount of flux recovered in the high resolution image versus
        the low resolution data. If the spectral axes don't match, one is
        interpolated onto the other.

        Parameters
        ----------
        plot : bool, optional
            Enable plotting.
        filename : str, optional
            Give filename for the plot save file. When specified, the plot
            is automatically saved.
        interp_to : 'lower' or 'upper', optional
            If the spectral axes don't match, interpolated onto the same.
            The default 'lower' interpolates to the spectral axis with the
            lowest resolution.
        '''

        # Add up the total intensity in the cubes and compare
        # Assumes that there is some masking such that noise
        # doesn't dominate

        if self.highres_convolved is not None:
            high_channel_intensity = \
                self.highres_convolved.sum(axis=(1, 2))
        else:
            high_channel_intensity = self.highres.sum(axis=(1, 2))
            Warning("Should run convolve_to_common before. Using unconvolved"
                    " cube.")

        if self.lowres_convolved is not None:
            low_channel_intensity = self.lowres_convolved.sum(axis=(1, 2))
        else:
            low_channel_intensity = self.lowres.sum(axis=(1, 2))
            Warning("Should run convolve_to_common before. Using unconvolved"
                    " cube.")

        # If the spectral axes are not the same, try interpolating to match
        if enable_interp:

            better_vres_high = \
                np.abs(self.lowres.header['CDELT3']) > \
                np.abs(self.highres.header['CDELT3'])

            # Invert which to interpolate to
            if interp_to == "higher":
                better_vres_high = not better_vres_high

            if better_vres_high:
                f = interp1d(np.round(self.highres.spectral_axis.value, 3),
                             high_channel_intensity.value)
                high_channel_intensity = f(np.round(self.lowres.spectral_axis.value, 3)) *\
                    self.highres.unit
            else:
                f = interp1d(np.round(self.lowres.spectral_axis.value, 3),
                             low_channel_intensity.value)
                low_channel_intensity = f(np.round(self.highres.spectral_axis.value, 3)) *\
                    self.lowres.unit


        self.fraction_flux_recovered = \
            high_channel_intensity.sum()/low_channel_intensity.sum()

        if plot:
            p.plot(self.highres.spectral_axis.value,
                   high_channel_intensity.value)

            p.plot(self.lowres.spectral_axis.value,
                   low_channel_intensity.value)

            p.xlabel("Spectral Axis ("+self.highres.spectral_axis.unit.to_string()+")")
            p.ylabel("Intensity ("+self.highres.unit.to_string()+")")

            if filename is None:
                p.show()
            else:
                p.savefig(filename)

        def run_all(self, highres_mask=None, lowres_mask=None, plot=True,
                    unit=u.K, freq=1420.40575177*u.MHz, filename=True,
                    use_dask=False, verbose=True):
            '''
            Run the complete comparisons.
            '''

            self.apply_mask(highres_mask=highres_mask, lowres_mask=lowres_mask)

            self.match_coords()

            self.convert_to(unit=unit, freq=freq)

            self.convolve_to_common(verbose=verbose, use_dask=use_dask)

            self.flux_recovered(plot=plot, filename=filename)


def _update_beam_in_hdr(hdr, beam):
    hdr["BMAJ"] = beam.major.value
    hdr["BMIN"] = beam.minor.value
    hdr["BPA"] = beam.pa.value
    return hdr


def padwithnans(vector, pad_width, iaxis, kwargs):
    vector[:pad_width[0]] = np.nan
    vector[-pad_width[1]:] = np.nan
    return vector


def auto_dask_map(cube, operation=convolve_fft, blocks=None, args=[],
                  kwargs={}, output_array=None, verbose=True):
    '''
    Based on the dimensions of the given cube, and the dimensions
    of the blocks, return an appropriate dask output.

    Currently, this only supports 2D and 3D inputs. For 2D blocks,
    an iterator is returned with 2D slices of the array. For 3D,
    a single dask array is returned.
    '''

    if len(blocks) == 2:

        if output_array is None:
            output_array = np.empty(cube.shape)

            def_slice = [slice(None)] * len(cube.shape)

        for chan, slice in dask_slice_iterator(cube, blocks):
            if verbose:
                print("On "+str(chan))
            def_slice[0] = slice(chan, chan+1)
            output_array[def_slice] = \
                slice.map_overlap(lambda a:
                                  operation(a, *args, **kwargs)).compute()

    else:
        dask_arr = da.from_array(cube.filled_data[:], blocks)

        if verbose:
            print("No channel iterations, so no print out.")

        # In the >2D case, map blocks to function in 2D planes

        output_array = \
            dask_arr.map_blocks(
                lambda a: _map_flatten(convolve_fft,
                                       a, args=args,
                                       kwargs=kwargs)).compute()

    return output_array


def dask_slice_iterator(cube, blocks):
    for chan in range(len(cube.spectral_axis)):
        yield chan, da.from_array(cube.filled_data[chan, :, :],
                                  blocks)

def _map_flatten(operation, arr, args=[], kwargs={}):

    restore_slice = [slice(None)] * len(arr.shape)

    for i in np.where(np.asarray(arr.shape)==1)[0]:
        restore_slice[i] = np.newaxis

    return operation(arr.squeeze(), *args, **kwargs)[restore_slice]
