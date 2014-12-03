
'''
Test some pruning
'''

from astrodendro import Dendrogram
from astropy.io.fits import getdata

img = getdata("hd22.13co.intintensity.fits")

d1 = Dendrogram.compute(img, verbose=True, min_delta=1.0, min_npix=2)

d2 = Dendrogram.compute(img, verbose=True, min_delta=0.5, min_npix=2)
d2.prune(min_delta=1.0)

from astrodendro.tests import test_pruning

test_pruning.compare_dendrograms(d1, d2)
