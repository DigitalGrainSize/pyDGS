# encoding: utf-8
"""
pyDGS - a Python framework for wavelet-based digital grain size analysis

pyDGS is an open-source project dedicated to provide a Python framework to
compute estimates of grain size distribution  using the continuous wavelet transform method
of Buscombe (2013) from an image of sediment where grains are clearly resolved.

This program implements the algorithm of:

Buscombe, D. (2013)
Transferable Wavelet Method for Grain-Size Distribution from Images of Sediment Surfaces and Thin Sections,
and Other Natural Granular Patterns. Sedimentology 60, 1709-1732

http://dbuscombe-usgs.github.io/docs/Buscombe2013_Sedimentology_sed12049.pdf

 Author:  Daniel Buscombe
           Marda Science, LLC
           Flagstaff, AZ
           daniel@mardascience.com
 First Revision January 18 2013

For more information visit https://github.com/dbuscombe-usgs/pyDGS
"""

# Written by Dr Daniel Buscombe, Marda Science LLC
#
# MIT License
#
# Copyright (c) 2020-22, Marda Science LLC
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import numpy as np
import sys, os
from imageio import imread
import pywt
from tqdm import tqdm
from skimage.restoration import denoise_wavelet, estimate_sigma
from functools import partial
# rescale_sigma=True required to silence deprecation warnings
_denoise_wavelet = partial(denoise_wavelet, rescale_sigma=True)
import scipy.stats as stats

# =========================================================
def rescale(dat,mn,mx):
    """
    rescales an input dat between mn and mx
    """
    m = min(dat.flatten())
    M = max(dat.flatten())
    return (mx-mn)*(dat-m)/(M-m)+mn

##====================================
def standardize(img):
    img = np.array(img)
    #standardization using adjusted standard deviation
    N = np.shape(img)[0] * np.shape(img)[1]
    s = np.maximum(np.std(img), 1.0/np.sqrt(N))
    m = np.mean(img)
    img = (img - m) / s
    img = rescale(img, 0, 1)
    del m, s, N

    return img

# =========================================================
# =========================================================
def dgs(image, resolution=1, maxscale=4, verbose=1, x=-0.5):

   if verbose==1:
      print("===========================================")
      print("======DIGITAL GRAIN SIZE: WAVELET==========")
      print("===========================================")
      print("=CALCULATE GRAIN SIZE-DISTRIBUTION FROM AN=")
      print("====IMAGE OF SEDIMENT/GRANULAR MATERIAL====")
      print("===========================================")
      print("======A PROGRAM BY DANIEL BUSCOMBE=========")
      print("====MARDASCIENCE, FLAGSTAFF, ARIZONA=======")
      print("========REVISION 4.2, APR 2022===========")
      print("===========================================")

   # exit program if no input folder given
   if not image:
      print('An image file is required!!!!!!')
      sys.exit(2)

   # print given arguments to screen and convert data type where necessary
   if image:
      print('Input image is '+image)
   #
   # if resolution:
   #    resolution = np.asarray(resolution,float)
   #    print('Resolution is '+str(resolution))
   #
   # if maxscale:
   #    maxscale = np.asarray(maxscale,int)
   #    print('Max scale as inverse fraction of data length: '+str(maxscale))
   #
   # if x:
   #    x = np.asarray(x, float)
   #    print('Area to volume conversion constant = '+str(x))

   # ======= stage 1 ==========================
   # read image
   if verbose==1:
      print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
      print('Processing image '+image)
   try:
       im = imread(image)   # read the image straight with imread
       im = np.squeeze(im)  # squeeze singleton dimensions
       if len(np.shape(im))>3:
           im = im[:, :, :3]            # only keep the first 3 bands

       if len(np.shape(im))==3: # if rgb, convert to grey
          im = (0.299 * im[:,:,0] + 0.5870*im[:,:,1] + 0.114*im[:,:,2]).astype('uint8')

       nx,ny = np.shape(im)
       if nx>ny:
          im=im.T

       im = standardize(im)

   except: # IOError:
       print('cannot open '+image)
       sys.exit(2)

   # # ======= stage 2 ==========================
   # Denoised image using default parameters of `denoise_wavelet`
   filter=False

   if filter:
      sigma_est = estimate_sigma(im, multichannel=False, average_sigmas=True)
      region = denoise_wavelet(im, multichannel=False, rescale_sigma=True,
                                 method='VisuShrink', mode='soft', sigma=sigma_est*2)
   else:
      region = im.copy()

   original = rescale(region,0,255)

   nx, ny = original.shape

   # ======= stage 3 ==========================
   # call cwt to get particle size distribution

   # P = []
   # for k in tqdm(np.linspace(1,nx-1,100)):
   #    [cfs, frequencies] = pywt.cwt(original[int(k),:], np.arange(1, ny/maxscale, 2),  'morl' , .5)
   #    period = 1. / frequencies
   #    power =(abs(cfs)) ** 2
   #    P.append(np.mean(np.abs(power), axis=1)/(period**2))

   # p = np.mean(np.vstack(P), axis=0)
   # p = np.array(p/np.sum(p))

   # #plt.plot(period, p,'m', lw=2); plt.show()

   # # get real scales by multiplying by resolution (mm/pixel)
   # scales = np.array(period)*resolution

   # ind = np.where(scales>2*np.pi)[0]
   # scales = scales[ind]
   # p = p[ind]
   # p = p/np.sum(p)

   P = []; M = []
   for k in np.linspace(1,nx-1,100):
      [cfs, frequencies] = pywt.cwt(original[int(k),:], np.arange(3, np.maximum(nx,ny)/maxscale, 1),  'morl' , .5)
      period = 1. / frequencies
      power =(abs(cfs)) ** 2
      power = np.mean(np.abs(power), axis=1)/(period**2)
      P.append(power)

      M.append(period[np.argmax(power)])

   p = np.mean(np.vstack(P), axis=0)
   p = np.array(p/np.sum(p))

   # get real scales by multiplying by resolution (mm/pixel)
   scales = np.array(period)*resolution

   srt = np.sqrt(np.sum(p*((scales-np.mean(M))**2)))

   # plt.plot(scales, p,'m', lw=2)

   p = p+stats.norm.pdf(scales, np.mean(M), srt/2)
   p = np.hstack([0,p])
   scales = np.hstack([0,scales])
   p = p/np.sum(p)

   # area-by-number to volume-by-number
   r_v = (p*scales**x) / np.sum(p*scales**x) #volume-by-weight proportion

   # ======= stage 5 ==========================
   # calc particle size stats

   pd = np.interp([.05,.1,.16,.25,.3,.5,.75,.84,.9,.95],np.hstack((0,np.cumsum(r_v))), np.hstack((0,scales)) )
   if verbose==1:
      print("d50 = "+str(pd[4]))

   mnsz = np.sum(r_v*scales)
   if verbose==1:
      print("mean size = "+str(mnsz))

   srt = np.sqrt(np.sum(r_v*((scales-mnsz)**2)))
   if verbose==1:
      print("stdev = "+str(srt))

   sk = (sum(r_v*((scales-mnsz)**3)))/(100*srt**3)
   if verbose==1:
      print("skewness = "+str(sk))

   kurt = (sum(r_v*((scales-mnsz)**4)))/(100*srt**4)
   if verbose==1:
      print("kurtosis = "+str(kurt))

   # plt.plot(scales, r_v,'k', lw=2); plt.show()

   # ======= stage 6 ==========================
   # return a dict object of stats
   return {'mean grain size': mnsz, 'grain size sorting': srt, 'grain size skewness': sk, 'grain size kurtosis': kurt, 'percentiles': [.05,.1,.16,.25,.3,.5,.75,.84,.9,.95], 'percentile_values': pd, 'grain size frequencies': r_v, 'grain size bins': scales}


# =========================================================
# =========================================================
if __name__ == '__main__':

   dgs(image, resolution=1, maxscale=8, verbose=0, x=-1)
