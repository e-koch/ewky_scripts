
import numpy as np
import skimage.morphology as mo
import skimage.measure as me
from skimage.segmentation import find_boundaries


def edge_curvature(mask, min_sep=5, average_over=3):
    '''
    Compute the menger curvature along the edges of the contours in the mask.
    '''

    labels = me.label(mask, neighbors=8, connectivity=2)

    edges = find_boundaries(labels, connectivity=2, mode='outer')

    pts = integer_boundaries(mask, edges, 0.5)

    curvature_mask = np.zeros_like(mask, dtype=float)

    for cont_pts in pts:
        # Last one is a duplicate
        cont_pts = cont_pts[:-1]

        num = cont_pts.shape[0]

        for i in xrange(num):

            curv = 0.0
            for j in xrange(min_sep, min_sep+average_over+1):
                curv += menger_curvature(cont_pts[i-j], cont_pts[i],
                                         cont_pts[(i+j) % num])

            y, x = cont_pts[i]

            if np.isnan(curv):
                curv = 0.0
            curvature_mask[y, x] = curv / average_over

    return curvature_mask


def menger_curvature(pt1, pt2, pt3, atol=1e-3):

    vec21 = np.array([pt1[0]-pt2[0], pt1[1]-pt2[1]])
    vec23 = np.array([pt3[0]-pt2[0], pt3[1]-pt2[1]])

    norm21 = np.linalg.norm(vec21)
    norm23 = np.linalg.norm(vec23)

    theta = np.arccos(np.dot(vec21, vec23)/(norm21*norm23))
    if np.isclose(theta-np.pi, 0.0, atol=atol):
        theta = 0.0

    dist13 = np.linalg.norm(vec21-vec23)

    return 2*np.sin(theta) / dist13


def integer_boundaries(mask, edges, level):
    '''
    Return the non-interpolated contour boundaries.
    '''

    all_pts = me.find_contours(mask, 0.5)

    int_pts = []

    for pts in all_pts:
        new_int_pts = np.zeros_like(pts, dtype=int)

        for i, pt in enumerate(pts):
            y, x = pt

            ceil = (np.ceil(y).astype(int), np.ceil(x).astype(int))
            floor = (np.floor(y).astype(int), np.floor(x).astype(int))

            if edges[ceil]:
                new_int_pts[i] = np.array(ceil)
            elif edges[floor]:
                new_int_pts[i] = np.array(floor)
            else:
                raise IndexError("Cannot find pixel in mask for " +
                                 str(pt))

        int_pts.append(new_int_pts)

    return int_pts


def curve(n, pts):
    '''
    The average curvature of the filament is found using the Menger curvature.
    The formula relates the area of the triangle created by the three points
    and the distance between the points. The formula is given as 4*area/|x-y||y-z||z-x|=curvature.
    The curvature is weighted by the Euclidean length of the three pixels.

    *Note:* The normalization is still an issue with this method. Its results
    should **NOT** be used.

    Parameters
    ----------

    n : int
        The number of the skeleton being analyzed.

    pts : list
          Contains the pixels contained in the inputted structure.

    Returns
    -------

    numer/denom : float
                  The value of the Menger Curvature.

    References
    ----------

    '''
    lenn = len(pts)
    kappa = []
    seg_len = []
    for i in range(lenn - 2):
        x1 = pts[i][0]
        y1 = pts[i][1]
        x2 = pts[i + 1][0]
        y2 = pts[i + 1][1]
        x3 = pts[i + 2][0]
        y3 = pts[i + 2][1]
        num = abs(2 * ((x2 - x1) * (y2 - y1) + (y3 - y2) * (x3 - x2)))
        den = np.sqrt((pow((x2 - x1), 2) + pow((y2 - y1), 2)) * (pow((x3 - x2),
                                                                     2) + pow((y3 - y2), 2)) * (pow((x1 - x3), 2) + pow((y1 - y3), 2)))
        if (den == 0):
            kappa.append(0)
        else:
            kappa.append(num / den)
        seg_len.append(
            fil_length(n, [[pts[i], pts[i + 1], pts[i + 2]]], initial=False)[0])
    numer = sum(kappa[i] * seg_len[i][0] for i in range(len(kappa)))
    denom = sum(seg_len[i][0] for i in range(len(seg_len)))
    if denom != 0:
        return numer / denom
    else:
        print n
        print pts
        raise ValueError('Sum of length segments is zero.')


def av_curvature(n, finalpix, ra_picks=100, seed=500):
    '''
    This function acts as a wrapper on curve. It calculates the average curvature
    by choosing 3 random points on the filament and calculating the curvature.
    The average of many iterations of this method is reported as the curvature
    for that skeleton.

    Parameters
    ----------

    n : int
        The number of the skeleton being analyzed.

    finalpix : list
               Contains the pixels contained in the inputted structure.

    ra_picks : int
               The number of iterations to run.

    seed : int
           Sets the seed.
    '''
    import numpy.random as ra
    seed = int(seed)
    ra.seed(seed=int(seed))
    ra_picks = int(ra_picks)

    curvature = []

    for i in range(len(finalpix)):
        if len(finalpix[i]) > 3:
            trials = []
            for _ in range(ra_picks):
                # REQUIRE NUMPY 1.7!!!
                picks = ra.choice(len(finalpix[i]), 3, replace=False)
                points = [finalpix[i][j] for j in picks]
                trials.append(curve(n, points))
            curvature.append(np.mean(trials))
        else:
            curvature.append("Fail")
    return curvature
