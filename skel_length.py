
'''

Length Testing

'''

import numpy as np
from skimage.morphology import label
import scipy.ndimage as nd
from copy import copy

struct2a = np.eye(2)
struct2b = struct2a[::-1]

struct3a = np.eye(3)
struct3b = struct3a[::-1]


def lengths(skeleton, verbose=False):
    '''
    Length finding via morphology.
    '''

    # 4-connected labels
    four_labels = label(skeleton, 4, background=0)

    four_sizes = nd.sum(skeleton, four_labels, range(np.max(four_labels) + 1))

    # Lengths is the number of pixels minus number of objects with more
    # than 1 pixel.
    four_length = np.sum(
        four_sizes[four_sizes > 1]) - len(four_sizes[four_sizes > 1])

    # Find pixels which a 4-connected and subtract them off the skeleton

    four_objects = np.where(four_sizes > 1)[0]

    skel_copy = copy(skeleton)
    for val in four_objects:
        skel_copy[np.where(four_labels == val)] = 0

    # Remaining pixels are only 8-connected
    # Lengths is same as before, multiplied by sqrt(2)

    eight_labels = label(skel_copy, 8, background=0)

    eight_sizes = nd.sum(
        skel_copy, eight_labels, range(np.max(eight_labels) + 1))

    eight_length = (
        (np.sum(eight_sizes) - 1) - np.max(eight_labels)) * np.sqrt(2)

    # If there are no 4-connected pixels, we don't need the hit-miss portion.
    if four_length == 0.0:
        conn_length = 0.0

    else:

        # Check 4 to 8-connected elements
        struct1 = np.array([[1, 0, 0],
                            [0, 1, 1],
                            [0, 0, 0]])

        struct2 = np.array([[0, 0, 1],
                            [1, 1, 0],
                            [0, 0, 0]])

        # Next check the three elements which will be double counted
        check1 = np.array([[1, 1, 0, 0],
                           [0, 0, 1, 1]])

        check2 = np.array([[0, 0, 1, 1],
                           [1, 1, 0, 0]])

        check3 = np.array([[1, 1, 0],
                           [0, 0, 1],
                           [0, 0, 1]])

        store = np.zeros(skeleton.shape)

        # Loop through the 4 rotations of the structuring elements
        for k in range(0, 4):
            hm1 = nd.binary_hit_or_miss(
                skeleton, structure1=np.rot90(struct1, k=k))
            store += hm1

            hm2 = nd.binary_hit_or_miss(
                skeleton, structure1=np.rot90(struct2, k=k))
            store += hm2

            hm_check3 = nd.binary_hit_or_miss(
                skeleton, structure1=np.rot90(check3, k=k))
            store -= hm_check3

            if k <= 1:
                hm_check1 = nd.binary_hit_or_miss(
                    skeleton, structure1=np.rot90(check1, k=k))
                store -= hm_check1

                hm_check2 = nd.binary_hit_or_miss(
                    skeleton, structure1=np.rot90(check2, k=k))
                store -= hm_check2

        conn_length = np.sqrt(
            2) * np.sum(np.sum(store, axis=1), axis=0)  # hits

    if verbose:
        print "Four Length: %s" % (four_length)
        print "Eight Length: %s" % (eight_length)
        print "Connect Length: %s" % (conn_length)

    return conn_length + eight_length + four_length


if __name__ == "__main__":

    import numpy as np
    import matplotlib.pyplot as p

    test_skel1 = np.zeros((10, 10))

    test_skel1[1, 1:3] = 1
    test_skel1[2, 3] = 1
    test_skel1[3, 4:6] = 1

    length1 = lengths(test_skel1)  # , verbose=True)
    print length1
    print 2 + np.sqrt(2) * 2
    print "Match: %s" % (length1 == 2 + np.sqrt(2) * 2)

    test_skel2 = np.eye(10)

    length2 = lengths(test_skel2)

    print length2
    print 9 * np.sqrt(2)
    print "Match: %s" % (length2 == 9 * np.sqrt(2))

    test_skel3 = np.zeros((10, 10))

    test_skel3[:, 5] = 1

    length3 = lengths(test_skel3)

    print length3
    print 9
    print "Match: %s" % (length3 == 9)

    test_skel4 = np.zeros((12, 12))

    test_skel4[0, 3] = 1
    test_skel4[1, 2] = 1
    test_skel4[2, 1] = 1
    test_skel4[3, 0] = 1
    test_skel4[4, 1] = 1
    test_skel4[5, 2] = 1
    test_skel4[6, 3] = 1
    test_skel4[7, 3] = 1
    test_skel4[8, 3] = 1
    test_skel4[9, 4] = 1
    test_skel4[10, 5] = 1
    test_skel4[11, 6] = 1
    test_skel4[11, 7] = 1
    test_skel4[11, 8] = 1
    test_skel4[10, 9] = 1
    test_skel4[9, 10] = 1
    test_skel4[9, 11] = 1

    length4 = lengths(test_skel4)  # , verbose=True)

    print length4
    print 11 * np.sqrt(2) + 5
    print "Match: %s" % (length4 == 11 * np.sqrt(2) + 5)

    test_skel5 = np.zeros((9, 12))

    test_skel5[8, 0] = 1
    test_skel5[8, 1] = 1
    test_skel5[7, 2] = 1
    test_skel5[7, 3] = 1
    test_skel5[6, 4] = 1
    test_skel5[5, 4] = 1
    test_skel5[4, 4] = 1
    test_skel5[3, 5] = 1
    test_skel5[2, 5] = 1
    test_skel5[1, 5] = 1
    test_skel5[0, 6] = 1
    test_skel5[0, 7] = 1
    test_skel5[1, 8] = 1
    test_skel5[1, 9] = 1
    test_skel5[2, 10] = 1
    test_skel5[3, 10] = 1
    test_skel5[4, 11] = 1
    test_skel5[5, 11] = 1

    length5 = lengths(test_skel5, verbose=True)

    print length5
    print 10 + 7 * np.sqrt(2)
    print "Match: %s" % (length5 == 10 + 7 * np.sqrt(2))

    test_skel6 = np.zeros((4, 4))

    test_skel6[1, 0:2] = 1
    test_skel6[2:4, 2] = 1

    length6 = lengths(test_skel6, verbose=True)

    print length6
    print 2 + np.sqrt(2)
    print "Match: %s" % (length6 == 2 + np.sqrt(2))
