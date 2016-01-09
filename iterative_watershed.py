
import skimage.morphology as mo
import skimage.measure as me
from skimage.feature import peak_local_max
import scipy.ndimage as nd
import numpy as np
import warnings

try:
    import cv2
    CV2_FLAG = True
except ImportError:
    warnings.warn("Cannot import cv2. Computing with scipy.ndimage")
    CV2_FLAG = False


def iterative_watershed(array, scale, start_value=5, end_value=3,
                        delta_value=1, mask_below=2):
    '''
    Iterative Watershed algorithm.
    '''

    initial_mask = array >= mask_below

    # initial_mask = remove_spurs(initial_mask, scale)

    initial_peaks = peak_local_max(array, min_distance=scale,
                                   threshold_abs=start_value)

    initial_markers = np.zeros_like(array, dtype=bool)
    initial_markers[initial_peaks[:, 0], initial_peaks[:, 1]] = True

    wshed_input = -array.copy()
    wshed_input[wshed_input > 0] = 0

    labels = mo.watershed(wshed_input, me.label(initial_markers),
                          mask=wshed_input < -2)

    # Now decrease the local maxima, trying to subdivide
    # regions that had a peak at the higher level.
    peak_levels = \
        np.arange(start_value-delta_value,
                  end_value-delta_value, -1*delta_value)
    for value in peak_levels:
        new_peaks = peak_local_max(array, min_distance=scale,
                                   threshold_abs=value)
        print(new_peaks)
        markers = initial_markers.copy()
        markers[new_peaks[:, 0], new_peaks[:, 1]] = True

        # Remove markers not in the last watershed
        markers *= labels > 0

        # Search for label regions that now have multiple peaks
        # and re-run the watershed on them
        for lab in range(1, labels.max()+1):
            num_peaks = np.sum(markers*(labels == lab))
            print(lab, num_peaks)
            if num_peaks == 1:
                continue
            elif num_peaks == 0:
                raise Exception("No peaks found??")
            else:
                split_marker = me.label(markers*(labels == lab))
                split_label = mo.watershed(wshed_input, split_marker,
                                           mask=labels == lab)
                for lab2 in range(2, split_label.max()+1):
                    labels[np.where(split_label == lab2)] = labels.max() + 1

    return labels


def remove_spurs(mask, min_distance=9):
    '''
    Remove spurious mask features with reconstruction.
    '''

    # Distance transform of the mask
    dist_trans = nd.distance_transform_edt(mask)

    # We don't want to return local maxima within the minimum distance
    # Use reconstruction to remove.
    seed = dist_trans + min_distance
    reconst = mo.reconstruction(seed, dist_trans, method='erosion') - \
        min_distance

    if CV2_FLAG:
        return cv2.morphologyEx((reconst > 0).astype("uint8"),
                                cv2.MORPH_DILATE,
                                mo.disk(min_distance).astype("uint8")).astype(bool)
    else:
        return mo.dilation(reconst > 0, selem=mo.disk(min_distance))
