
from spectral_cube import SpectralCube
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
import matplotlib.pyplot as p

# hi_cube = SpectralCube.read(
#     '/media/eric/MyRAID/M33/14B-088/HI/imaging/south_arm_800_1200.image.fits')

hi_cube = SpectralCube.read(
    '/media/eric/MyRAID/M33/14B-088/HI/full_imaging/M33_14B-088_HI.clean.image.fits',
    mode='denywrite')

# hi_cube = hi_cube.to(u.K, hi_cube.beam.jtok_equiv(1420.40575177*u.MHz))

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

posn = SkyCoord("01h33m21.287", "+30d32m16.110", frame='icrs')


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
hi_spectrum = hi_spectrum.to(u.K, hi_cube.beam.jtok_equiv(1420.40575177*u.MHz))

co21_spectrum = co21_cube[:, co21_posn[0], co21_posn[1]]

plot_posns = False

if plot_posns:
    ax = p.subplot(131)
else:
    p.figure(figsize=(8, 4))
    ax = p.subplot(111)
    ax.plot(hi_spectrum.spectral_axis/1000., hi_spectrum.value,
            drawstyle='steps')
    ax.set_xlabel('Velocity (km/s)')
    ax.set_ylabel('HI Surface Brightness (K)')
    for tick in ax.get_xticklabels():
        tick.set_rotation(45)

    # Fill in HISA regions b/w spectra

    co21_min = co21_cube.spectral_axis[10]
    co21_mid = co21_cube.spectral_axis[9]
    co21_max = co21_cube.spectral_axis[8]

    hi_posn_min = hi_cube.closest_spectral_channel(co21_min)+1
    hi_posn_mid = hi_cube.closest_spectral_channel(co21_mid)
    hi_posn_max = hi_cube.closest_spectral_channel(co21_max)-1

    # co21_fillbetween = np.empty((hi_posn_min-hi_posn_max,))

    # co21_fillbetween[:hi_posn_mid-hi_posn_max] = co21_spectrum.value[9]
    # co21_fillbetween[hi_posn_mid-hi_posn_max:] = co21_spectrum.value[10]

    # y1 = hi_spectrum.value[hi_posn_max:hi_posn_min]
    # y2 = co21_fillbetween * 525  # scale between ax and ax_2

    # ax.fill_between(hi_spectrum.spectral_axis[hi_posn_max:hi_posn_min]/1000.,
    #                 y1, y2, where=y2 >= y1, facecolor='k', interpolate=False,
    #                 alpha=0.25)

    y1 = np.array([-10]*27)
    y2 = np.array([50]*27)

    ax.fill_between(hi_spectrum.spectral_axis[hi_posn_max:hi_posn_min]/1000.,
                    y1, y2, where=y2 >= y1, facecolor='k', alpha=0.25)

    ax_2 = ax.twinx()
    ax_2.plot(co21_spectrum.spectral_axis/1000., co21_spectrum.value, 'g',
              drawstyle='steps')
    ax_2.set_ylabel('CO(2-1) Surface Brightness (K)')

if plot_posns:
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

p.tight_layout()

# p.savefig("hi_co_overlay_"+new_coord.to_string('hmsdms').replace(' ', '_')+".pdf")
