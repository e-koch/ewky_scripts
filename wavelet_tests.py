
'''
Ensure that the turbustat form gives the same slope as using the astropy
kernel (modulo some normalization factors...)
'''

from turbustat.statistics.wavelets.wavelet_transform import Mexican_hat

from astropy.convolution import MexicanHat2DKernel

import numpy as np
import matplotlib.pyplot as p

# 256x256 grid
y = np.arange(-128, 128)
x = np.arange(-128, 128)

k = np.fft.fftfreq(256)
l = np.fft.fftfreq(256)

width = 10.

turb_hat = Mexican_hat()
psi = turb_hat.psi(y/width, x/width) / (2*np.pi*width**4.)
psi_ft = turb_hat.psi_ft(2*np.pi*width*k, 2*np.pi*width*l) / (2*np.pi*width**4.)

apy_hat = MexicanHat2DKernel(10, x_size=256, y_size=256)

# p.ioff()

# p.subplot(131)
# p.title("Astropy")
# p.imshow(apy_hat)
# p.subplot(132)
# p.title("Turb")
# p.imshow(psi)
# p.subplot(133)
# p.title("Difference")
# p.imshow(apy_hat.array - psi)
# p.show()

p.subplot(131)
p.title("Astropy")
p.imshow(apy_hat)
p.subplot(132)
p.title("Turb FT")
p.imshow(np.fft.fftshift(np.fft.fft2(psi_ft)).real)
p.subplot(133)
p.title("Difference")
p.imshow(apy_hat.array - np.fft.fftshift(np.fft.fft2(psi_ft)).real)
p.show()

# p.ion()