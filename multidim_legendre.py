#!/usr/bin/python

'''
2D Legendre Fitting

Based on Numpy legendre routines

'''

import numpy as np
from numpy.polynomial.legendre import legfit,legvander,legvander2d
from numpy.linalg import lstsq
import matplotlib.pyplot as p
import warnings

def leg2dfit(array, deg, rcond=None, full=False, w=None):
  '''

  2D Legendre Fitting.
  This routine has an issue with mapping back off of the [-1,1] interval. The issue lies in the normalization.

  '''

  order = int(deg) + 1
  flat_array = array.ravel()

  # find positions of nans
  goodpts = np.where(np.isnan(flat_array))

  if w is not None:
    if w.ndim != 2:
        raise TypeError("expected 2D array for w")
    if array.shape != w.shape:
        raise TypeError("expected array and w to have same shape")
    goodpts = np.where(~np.isnan(flat_array+w.ravel()))
    w = w.ravel()[goodpts]

  x, y = np.meshgrid(np.linspace(-1,1,array.shape[0]), np.linspace(-1,1,array.shape[1]))
  x, y = x.ravel()[goodpts], y.ravel()[goodpts]

  # check arguments.
  if deg < 0 :
      raise ValueError("expected deg >= 0")

  # set up the least squares matrices in transposed form
  lhs = legvander2d(x, y, [deg,deg]).T
  rhs = flat_array[goodpts].T

  if w is not None:
    lhs = lhs * w
    rhs = rhs * w

  # set rcond
  if rcond is None :
    rcond = len(x)*np.finfo(x.dtype).eps

  # scale the design matrix and solve the least squares equation
  if issubclass(lhs.dtype.type, np.complexfloating):
      scl = np.sqrt((np.square(lhs.real) + np.square(lhs.imag)).sum(1))
  else:
      scl = np.sqrt(np.square(lhs).sum(1))
      print scl

  c, resids, rank, s = lstsq(lhs.T/scl, rhs.T, rcond)
  c = (c.T/scl).T

  # warn on rank reduction
  if rank != order**2. and not full:
      msg = "The fit may be poorly conditioned"
      warnings.warn(msg)
  if full :
    return c, [resids, rank, s, rcond]
  else :
    return c

def full_leg2dfit(array,deg,patch_size,pix_to_pc,full=True):
  if len(deg)!=2: raise ValueError("deg must be in form [xdeg,ydeg]")
  if patch_size%2==0: raise ValueError("patch_size must be odd number")
  size_orig = array.shape
  patch_offset = int((patch_size-1)/2.)
  # Pad array with nans so edges can be properly fit (routine already ignores nans)
  arr = pad(array,patch_offset,padwithnans)

  # Add edges of array to padded regions
  arr[:patch_offset,patch_offset:arr.shape[1]-patch_offset] = array[:patch_offset,:]
  arr[arr.shape[0]-patch_offset:,patch_offset:arr.shape[1]-patch_offset] = array[array.shape[0]-patch_offset:,:]
  arr[patch_offset:arr.shape[0]-patch_offset,:patch_offset] = array[:,:patch_offset]
  arr[patch_offset:arr.shape[0]-patch_offset,arr.shape[1]-patch_offset:] = array[:,array.shape[1]-patch_offset:]

  size_new = arr.shape
  # Set up resulting arrays
  grad_x = np.zeros((size_orig))
  # grad_x_err = np.zeros((size_orig))
  grad_y = np.zeros((size_orig))
  grad_x_2 = np.zeros((size_orig))
  grad_y_2 = np.zeros((size_orig))

  len_x = patch_size
  len_y = patch_size
  x = np.linspace(-1,1,len_x)
  y = np.linspace(-1,1,len_y)
  meshx,meshy = np.meshgrid(x,y)

  for i in range(patch_offset+5,size_new[0]-(patch_offset+5)):
    for j in range(patch_offset+5,size_new[1]-(patch_offset+5)):
      missing = np.isnan(arr[i-patch_offset:i+patch_offset,j-patch_offset:j+patch_offset].ravel())
      mask_arr = arr[i-patch_offset:i+patch_offset,j-patch_offset:j+patch_offset].ravel()[~missing]

      mask_meshx = meshx.ravel()[~missing]
      mask_meshy = meshy.ravel()[~missing]
      if len(mask_meshy)==0: #Has to hold for both ### For catching negligible regions before fitting routine
        fit = [np.empty((3,3)),np.NaN]
        fit[0][:,:] = np.NaN
        # fit[1][0][:,:] = np.NaN
        coef = np.reshape(fit[0],(deg[0]+1,deg[1]+1))
        # resid = np.reshape(fit[0],(deg[0]+1,deg[1]+1))
        #print "No points to fit :(%s,%s)" % (i,j)

      else:
        try:
          #raise np.linalg.linalg.LinAlgError
          fit = leg2dfit(mask_meshx,mask_meshy,mask_arr,deg,full=True)
        except np.linalg.linalg.LinAlgError:
          print "Failed to fit :(%s,%s)" % (i,j)
          fit = [np.empty((3,3)),np.NaN]
          fit[0][:,:] = np.NaN
          print fit[1][0]
        coef = np.reshape(fit[0],(deg[0]+1,deg[1]+1))
        # resid = np.reshape(fit[1][0],(deg[0]+1,deg[1]+1))
        #print fit[0],coef
        if coef[0,0]>1.2*arr[i,j] or coef[0,0]<0.8*arr[i,j]:
          coef = coef#np.tril(coef[::-1])[::-1]
          #print "Fit fairly off: (%s,%s)" % (i,j)
          # grad_x[i-patch_offset,j-patch_offset] = np.NaN
          # grad_y[i-patch_offset,j-patch_offset] = np.NaN
          # grad_x_2[i-patch_offset,j-patch_offset] = np.NaN
          # grad_y_2[i-patch_offset,j-patch_offset] = np.NaN
        else:
          coef = coef#np.tril(coef[::-1])[::-1]
      #fit_arr = leggrid2d(x,y,coef).T * ~np.isnan(arr)

      #Input fits into Result arrays
        grad_x[i-patch_offset,j-patch_offset] = coef[0,1] * (1/float(pix_to_pc)) * patch_offset**-1
        # grad_x_err[i-patch_offset,j-patch_offset] = resid[0,1] * (1/float(pix_to_pc)**2.) * patch_offset**-2
        grad_y[i-patch_offset,j-patch_offset] = coef[1,2] * (1/float(pix_to_pc)) * patch_offset**-1
        grad_x_2[i-patch_offset,j-patch_offset] = coef[0,2] #* pix_to_pc * patch_offset
        grad_y_2[i-patch_offset,j-patch_offset] = coef[2,2] #* pix_to_pc * patch_offset

  grad_mag = np.sqrt(grad_x**2. + grad_y**2.)
  grad_dirn = np.arctan2(grad_y,grad_x)

  # Return the fit array, coefficients of the fit, and the residuals etc. from lstsq
  return grad_x, grad_y, grad_mag, grad_dirn, grad_x_2, grad_y_2

def legfunc_transform(coef,plot=False):
  # Given an array of legendre coefficients, plots the sum of
  # each degree of coefficients
  xdeg = coef.shape[0]-1
  ydeg = coef.shape[1]-1
  if xdeg!=ydeg: raise TypeError("X and Y degrees must be the same.")
  coef_sq = abs(coef)
  #coef_sq = coef_sq[::-1] #Flips the array, sum can then be found from the trace
  norm = range(1,xdeg+2)# + range(1,xdeg+1)[::-1]
  sum_coef = []
  for i in range(-xdeg,1):#xdeg +1):
    sum_coef.append(np.trace(coef_sq[::-1],offset=i))
  sum_coef = [sum_coef[i]/norm[i] for i in range(len(sum_coef))]
  coef_range = range(xdeg+1)
  if plot:
    p.plot(coef_range,sum_coef)
    p.show()
  else: pass
