#!/usr/bin/python

import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as p
from matplotlib import rc

rc('text', usetex=True)


def func1(x, k, t):
    term1 = 1. / ((k - 4.) ** 2. + 1.)
    term2 = np.cos(k * x + k ** 2. / 4. - 10. * k ** (0.8) * t)
    return term1 * term2


def func2(x, k, t):
    term1 = 1. / ((k - 4.) ** 2. + 1.)
    term2 = np.sin(k * x + k ** 2. / 4. - 10. * k ** (0.8) * t)
    return term1 * term2

x = np.linspace(-10, 10, 2001)

psi1 = []
psi2 = []
for i in x:
    val = quad(lambda k: func1(i, k, 0), 3., 5.)[
        0] ** 2. + quad(lambda k: func2(i, k, 0), 3., 5.)[0] ** 2.
    if val >= 0.49 and val <= 0.51:
        print i
    psi1.append(val / 2.5)
    val2 = quad(lambda k: func1(i, k, 1), 3.5, 4.5)[
        0] ** 2. + quad(lambda k: func2(i, k, 1), 3.5, 4.5)[0] ** 2.
    if val2 >= 0.17 and val2 <= 0.19:
        print i
    psi2.append(val2 / 2.5)


p.plot(x, psi1, "k", x, psi2, "b", label="test")
p.grid(True)
p.title("Wavepacket at t-=0 and t=1")
p.xlabel("Time")
p.ylabel("$|\Psi|$")
p.show()
