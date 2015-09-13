
from astropy.analytic_functions import blackbody_nu
from astropy import units as u
import numpy as np

'''
Dust continuum to mass.

Given an intensity, calculates the dust mass assuming constant temperature,
opacity, etc, etc...
'''

nu = 230.538 * u.GHz

T_dust = 15 * u.K
kappa_dust = 0.6 * u.cm**2/u.g

# Assume the Canonical value here
dust_to_gas = 100.

beamsize = 1 * u.arcsec

# Assuming circular here. Drop units to avoid 'beam' units in intensity
beamarea = (np.pi * beamsize ** 2).to(u.sr)

# In mJy/bm, but astropy isn't decomposing automatically
intensity = (253.54 * 1e-3) * 1e-26 * u.erg/(u.s * u.Hz * u.cm**2)

intensity = intensity / beamarea

Sigma_gas = (dust_to_gas / kappa_dust) * \
    (intensity / blackbody_nu(nu, T_dust))

Sigma_gas = Sigma_gas.to(u.Msun/u.pc**2)

# Multiply by physical beamwidth size
# In the case of M33, 1" = 4 pc
phys_beamsize = (4 * u.pc)
phys_beamarea = (np.pi*phys_beamsize**2)

dustmass = Sigma_gas * phys_beamarea

print "Mass Surface Density: %s" % (Sigma_gas)
print "Mass in beam: %s" % (dustmass)
