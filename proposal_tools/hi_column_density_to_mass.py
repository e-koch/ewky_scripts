
import astropy.units as u
import astropy.constants as con
import numpy as np
from radio_beam import Beam

'''
Recording the conversion from NHI to M, since I keep getting close to the
apparent correct numerical factor with slight differences everytime I try.

In the optically-thin limit, NHI = 1.82e18 * int T_B dv, where T_B is in K and
dv is in km/s. The factor of 1.82e18 must then have units of cm^-2 / K km s^-1.

To convert to mass, M = m_H * int NHI dA = m_H * NHI * D^2 * theta^2 * pi/4.

The relation in terms of the flux density is M = 2.36e5 * D^2 * int S_v dv
where 2.36e5 Msol / Mpc^2 * Jy * km s^-1.

There aren't many places where I've found the direct conversion of mass from
the line brightness (K km s^-1). Two editions of Tools of Radio Astro have
two different values: 5.30e2 and 0.39e2. Neither of these can be correct,
since the scaling from Jy to K is ~1e5.
This equation is M = E * D^2 theta^2 * int T_B dv, with theta in arcsec.

Here I recover the typically quoted 2.36e5 conversion factor, and using that,
get 0.34 Msol / Mpc^2 * K * km s^-1 * arcsec^2.

'''

nu = 1.420405 * u.GHz

C = 1.8224e18 * (u.cm**-2 / (u.K * u.km * u.s**-1))

D = C * con.m_p.to(u.Msun) * \
    (Beam(1 * u.rad).jtok(nu) / u.Jy) * (1 * u.cm / (1 * u.cm).to(u.Mpc))**2

# Now the beam is assumed to NOT be given the FWHM parameters. So I think we
# need to convert out 1 factor of the FWHM on the area
fwhm_area_factor = np.pi / (4 * np.log(2))

D *= fwhm_area_factor

true_D = 2.36e5 * D.unit

frac_diff = 100 * np.abs(D.value - true_D.value) / true_D.value

print("Fraction difference: {}".format(frac_diff))
# Off by ~0.5%! Success!

# Now what's the conversion factor if we don't go back to Jy and stay in K.
# To match the usual scaling, use a 1 arcsecond beam
E = D / ((Beam(1 * u.arcsec).jtok(nu) / u.Jy) * fwhm_area_factor) / u.arcsec**2

print("Using flux density: {0} {1}".format(D.value, D.unit.to_string()))
print("Using brightness temp: {0} {1}".format(E.value, E.unit.to_string()))
