
'''
Wrap a given function into a pool
'''

from multiprocessing import Pool
from itertools import izip, repeat


def pool_process(func, args=None, ncores=1):

    def single_input(a):
        return func(*a)

    # Check for inputs which need to be repeated.
    for arg in args:
        try:
            if len(arg) == 1:
                repeat_it = True
        except TypeError:
            repeat_it = True

        if repeat_it:
            arg = repeat(arg)

    pool = Pool(pool_process=ncores)
    output = pool.map(single_input, izip(*args))

    pool.close()
    pool.join()

    return output

