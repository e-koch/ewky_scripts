
from astropy import units as u
from astropy.constants import c
import numpy as np


def optical_to_radio(vel, freq, rest):

    return freq * vel / rest


def radio_to_optical(vel, freq, rest):

    return rest * vel / freq


def test_opt_to_rad():

    opt_vel = -2126.1876453900204 * (u.km/u.s)

    rad_vel = -2141.374699999949 * (u.km/u.s)

    rest = 1.4 * u.GHz

    freq = 1.41 * u.GHz

    np.testing.assert_almost_equal(opt_vel.value,
                                   radio_to_optical(rad_vel, freq, rest).value)

    np.testing.assert_almost_equal(rad_vel.value,
                                   optical_to_radio(opt_vel, freq, rest).value)


def direct_radio_to_optical(velocity):
    try:
        velocity.units
        unit_flag = True
    except AttributeError:
        unit_flag = False

    if unit_flag:
        return (c * velocity) / (c - velocity)

    else:
        return (c.value * velocity) / (c.value - velocity)


def direct_optical_to_radio(velocity):
    try:
        velocity.units
        unit_flag = True
    except AttributeError:
        unit_flag = False

    if unit_flag:
        return (c * velocity) / (c + velocity)

    else:
        return (c.value * velocity) / (c.value + velocity)

