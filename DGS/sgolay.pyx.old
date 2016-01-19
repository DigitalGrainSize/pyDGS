"""
 Author:  Daniel Buscombe
           Grand Canyon Monitoring and Research Center
           United States Geological Survey
           Flagstaff, AZ 86001
           dbuscombe@usgs.gov
 First Revision January 18 2013

For more information visit https://github.com/dbuscombe-usgs/DGS-python

    GNU Lesser General Public License, Version 3
    (http://www.gnu.org/copyleft/lesser.html)
    
    This software is in the public domain because it contains materials that
    originally came from the United States Geological Survey, an agency of the
    United States Department of Interior. For more information, 
    see the official USGS copyright policy at
    http://www.usgs.gov/visual-id/credit_usgs.html#copyright
    Any use of trade, product, or firm names is for descriptive purposes only 
    and does not imply endorsement by the U.S. government.
"""

from __future__ import division
import numpy as np
cimport numpy as np
cimport cython
import scipy.signal as sp

# =========================================================
cdef class sgolay2d:
            
    cdef object data
    
    # =========================================================
    @cython.boundscheck(False)
    @cython.cdivision(True)
    @cython.wraparound(False)
    @cython.nonecheck(False)
    def __init__(self, np.ndarray[np.uint8_t, ndim=2] z, int window_size, int order ):
       """
       do 2d filtering on matrix
       from http://www.scipy.org/Cookbook/SavitzkyGolay
       """

       # number of terms in the polynomial expression
       cdef float n_terms = ( order + 1 ) * ( order + 2)  / 2.0

       if  window_size % 2 == 0:
           raise ValueError('window_size must be odd')

       if window_size**2 < n_terms:
           raise ValueError('order is too high for the window size')

       cdef int half_size = window_size // 2

       # exponents of the polynomial. 
       # p(x,y) = a0 + a1*x + a2*y + a3*x^2 + a4*y^2 + a5*x*y + ... 
       # this line gives a list of two item tuple. Each tuple contains 
       # the exponents of the k-th term. First element of tuple is for x
       # second element for y.
       # Ex. exps = [(0,0), (1,0), (0,1), (2,0), (1,1), (0,2), ...]
       cdef list exps = [] 
       exps = [ (k-n, n) for k in xrange(order+1) for n in xrange(k+1) ]

       # coordinates of points
       cdef np.ndarray[np.float64_t,ndim=1] ind = np.arange(-half_size, half_size+1, dtype=np.float64)
       cdef np.ndarray[np.float64_t,ndim=1] dx = np.empty(len(ind)*window_size, dtype=np.float64) 
       dx = np.repeat( ind, window_size )
       cdef np.ndarray[np.float64_t,ndim=1] dy = np.empty(len(ind)*window_size, dtype=np.float64) 
       dy = np.tile( ind, [window_size, 1]).reshape(window_size**2, )

       # build matrix of system of equation
       cdef np.ndarray[np.float64_t,ndim=2] A = np.empty( (window_size**2, len(exps)), dtype=np.float64 )
       for i, exp in enumerate( exps ):
          A[:,i] = (dx**exp[0]) * (dy**exp[1])

       # pad input array with appropriate values at the four borders
       cdef np.ndarray[np.float64_t,ndim=2] Z = np.zeros( (z.shape[0] + 2*half_size, z.shape[1] + 2*half_size), dtype=np.float64 )
       # top band
       Z[:half_size, half_size:-half_size] =  z[0, :] -  np.abs( np.flipud( z[1:half_size+1, :] ) - z[0, :] )
       # bottom band
       Z[-half_size:, half_size:-half_size] = z[-1, :]  + np.abs( np.flipud( z[-half_size-1:-1, :] )  -z[-1, :] )
       # left band
       Z[half_size:-half_size, :half_size] = np.tile( z[:,0].reshape(-1,1), [1,half_size]) - np.abs( np.fliplr( z[:, 1:half_size+1] ) - np.tile( z[:,0].reshape(-1,1), [1,half_size]) )
       # right band
       Z[half_size:-half_size, -half_size:] =  np.tile( z[:,-1].reshape(-1,1), [1,half_size] ) + np.abs( np.fliplr( z[:, -half_size-1:-1] ) - np.tile( z[:,-1].reshape(-1,1), [1,half_size] ) )
       # central band
       Z[half_size:-half_size, half_size:-half_size] = z

       # top left corner
       Z[:half_size,:half_size] = z[0,0] - np.abs( np.flipud(np.fliplr(z[1:half_size+1,1:half_size+1]) ) - z[0,0] )
       # bottom right corner
       Z[-half_size:,-half_size:] = z[-1,-1] + np.abs( np.flipud(np.fliplr(z[-half_size-1:-1,-half_size-1:-1]) ) - z[-1,-1] )

       # top right corner
       Z[:half_size,-half_size:] = Z[half_size,-half_size:] - np.abs( np.flipud(Z[half_size+1:2*half_size+1,-half_size:]) - Z[half_size,-half_size:] )
       # bottom left corner
       Z[-half_size:,:half_size] = Z[-half_size:,half_size].reshape(-1,1) - np.abs( np.fliplr(Z[-half_size:, half_size+1:2*half_size+1]) - Z[-half_size:,half_size].reshape(-1,1) )

       cdef np.ndarray[np.float64_t,ndim=2] m = np.zeros((window_size, window_size), dtype=np.float64)
       # solve system and convolve
       m = np.linalg.pinv(A)[0].reshape((window_size, -1))
       cdef np.ndarray[np.float64_t,ndim=2] out = np.zeros( np.shape(z), dtype=np.float64 )
       out = sp.fftconvolve(Z.astype('f'), m.astype('f'), mode='valid')
       self.data = out
       return

    # =========================================================    
    def getdata(self):
       return self.data
       
       

