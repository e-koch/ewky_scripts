
'''
Extract PV slices of the 'arm' coming of NGC-1333.
'''

from pvextractor import extract_pv_slice
from pvextractor.geometry import Path
from spectral_cube import SpectralCube, BooleanArrayMask
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import matplotlib.pyplot as p

# Start with the full COMPLETE Perseus cube
# Missing end card in the header, so we have to manually open it

hdu = fits.open("/srv/astro/erickoch/enzo_sims/complete/PerA_13coFCRAO_F_xyv.fits",
                ignore_missing_end=True)

mask = BooleanArrayMask(mask=np.isfinite(hdu[0].data), wcs=WCS(hdu[0].header))

sc = SpectralCube(data=hdu[0].data, wcs=WCS(hdu[0].header),
                  mask=mask)

# Select the region of interest
sc_small = sc[120:300, 156:346, 677:958]

# Ends for the PV slice
# ends = [(107, 74), (151, 76), (220, 54)]
ends = [(105, 70), (145, 80), (227, 28), (238, 1)]

xy = Path(ends, width=10)

pv = extract_pv_slice(sc_small, xy)

p.subplot(121)
p.imshow(sc_small.moment0().value, cmap="binary", origin="lower")
p.plot([i for i, j in ends], [j for i, j in ends], 'b-')

p.subplot(122)
p.imshow(pv.data, cmap="binary", origin="lower")
p.colorbar()

p.show()

hdu.close()

# Now examine a subsection of the region that overlaps with the edge of
# NGC-1333. I'm not sure the instrument that took the data, but it maps
# 13CO(1-0) too.

sc = SpectralCube.read("/srv/astro/surveys/classy/N1333.13co.fits")
sc = sc[150:300, :, :]

ends = [(57, 72), (73, 83), (105, 85)]

xy = Path(ends, width=10)

pv = extract_pv_slice(sc, xy)

p.subplot(121)
p.imshow(sc.moment0().value, cmap="binary", origin="lower")
p.plot([i for i, j in ends], [j for i, j in ends], 'b-')

p.subplot(122)
p.imshow(pv.data, cmap="binary", origin="lower")
p.colorbar()

p.show()
