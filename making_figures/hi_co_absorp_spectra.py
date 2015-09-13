
from spectral_cube import SpectralCube
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
import matplotlib.pyplot as p

hi_cube = SpectralCube.read(
    '/media/eric/MyRAID/M33/14B-088/HI/imaging/south_arm_800_1200.image.fits')

hi_cube = hi_cube.to(u.K, hi_cube.beam.jtok_equiv(1420.40575177*u.MHz))

co21_cube = SpectralCube.read('/media/eric/Data_3/M33/IRAM/m33.co21_iram.fits')
co21_cube = co21_cube.spectral_slab(*hi_cube.spectral_extrema)

# Only keep the overlap
lat_extrema = co21_cube.latitude_extrema
long_extrema = co21_cube.longitude_extrema
hi_cube = hi_cube.subcube(xlo=long_extrema[1],
                          xhi=long_extrema[0],
                          ylo=lat_extrema[0],
                          yhi=lat_extrema[1])

# Now find the position in each cube corresponding to a given
# coordinate

posn = SkyCoord("1h33m25.684", "+30d32m17.009", frame='icrs')


def get_closest_posn(posn, spatial_footprint):
    '''
    '''

    spatial_coords = SkyCoord(ra=spatial_footprint[1],
                              dec=spatial_footprint[0],
                              frame='icrs')

    min_posn = spatial_coords.separation(posn).argmin()

    twod_posn = np.unravel_index(min_posn, spatial_footprint[0].shape)

    return twod_posn

hi_posn = get_closest_posn(posn, hi_cube.spatial_coordinate_map)

co21_posn = get_closest_posn(posn, co21_cube.spatial_coordinate_map)

hi_spectrum = hi_cube[:, hi_posn[0], hi_posn[1]]

co21_spectrum = co21_cube[:, co21_posn[0], co21_posn[1]]


p.plot(hi_spectrum.spectral_axis, hi_spectrum.value)
p.plot(co21_spectrum.spectral_axis, co21_spectrum.value)