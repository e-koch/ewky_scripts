
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

posn = SkyCoord("1h33m14.905", "+30d33m17.479", frame='icrs')


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

ax = p.subplot(131)
ax.plot(hi_spectrum.spectral_axis/1000., hi_spectrum.value)
ax.set_xlabel('Velocity (km/s)')
ax.set_ylabel('HI Surface Brightness (K)')
for tick in ax.get_xticklabels():
    tick.set_rotation(45)
ax_2 = ax.twinx()
ax_2.plot(co21_spectrum.spectral_axis/1000., co21_spectrum.value, 'g')
ax_2.set_ylabel('CO(2-1) Surface Brightness (K)')

p.subplot(233)
hi_mom0 = np.nansum(hi_cube.filled_data[:].value, axis=0)
hi_slice = [slice(hi_posn[0]-30, hi_posn[0]+30),
            slice(hi_posn[1]-30, hi_posn[1]+30)]
p.imshow(hi_mom0[hi_slice], origin='lower')
p.scatter([30], [30], c='r', s=40)

p.subplot(236)
co21_mom0 = co21_cube.moment0().value
co21_slice = [slice(co21_posn[0]-30, co21_posn[0]+30),
              slice(co21_posn[1]-30, co21_posn[1]+30)]
p.imshow(co21_mom0[co21_slice], origin='lower')
p.scatter([30], [30], c='r', s=40)

p.savefig("hi_co_overlay_"+new_coord.to_string('hmsdms').replace(' ', '_')+".pdf")
