
'''
OpenCV is nice and optimized. ndimage is a little less so.
I just want to test if the outputs match, and the data types needed for that.
'''

import cv2
from scipy import ndimage as nd
import numpy as np
from skimage.morphology import disk


def tophat_test(plot=True):

    arr = np.ones((100, 100), dtype="float64")

    yy, xx = np.mgrid[-50:50, -50:50]

    arr[xx**2 + yy**2 <= 30**2] = 0.0

    struct = disk(32).astype("uint8")

    bth_nd = nd.black_tophat(arr, structure=struct)

    bth_cv = cv2.morphologyEx(arr, cv2.MORPH_BLACKHAT, struct)

    if plot:
        import matplotlib.pyplot as p
        p.subplot(131)
        p.imshow(bth_nd)
        p.colorbar()
        p.subplot(132)
        p.imshow(bth_cv)
        p.colorbar()
        p.subplot(133)
        p.imshow(np.abs(bth_nd-bth_cv))
        p.colorbar()
        p.show()

    np.testing.assert_equal(bth_nd, bth_cv)


if __name__ == "__main__":

    tophat_test()
