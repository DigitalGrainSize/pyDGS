# encoding: utf-8
"""
pyDGS - a Python framework for wavelet-based digital grain size analysis

pyDGS is an open-source project dedicated to provide a Python framework to compute estimates of grain size distribution  using the continuous wavelet transform method of Buscombe (2013) from an image of sediment where grains are clearly resolved. DOES NOT REQUIRE CALIBRATION

This program implements the algorithm of:

Buscombe, D. (2013) Transferable Wavelet Method for Grain-Size Distribution from Images of Sediment Surfaces and Thin Sections, and Other Natural Granular Patterns. Sedimentology 60, 1709-1732

http://dbuscombe-usgs.github.io/docs/Buscombe2013_Sedimentology_sed12049.pdf

 Author:  Daniel Buscombe
           Grand Canyon Monitoring and Research Center
           United States Geological Survey
           Flagstaff, AZ 86001
           dbuscombe@usgs.gov
 Revision Mar 1, 2016
 First Revision January 18 2013   

For more information visit https://github.com/dbuscombe-usgs/pyDGS

:install:
    python setup.py install
    sudo python setup.py install
    pip install pyDGS
    
:test:
    python -c "import DGS; DGS.test.dotest_web()"

    python
    import DGS

    image_file = '/home/sed_images/my_image.png'

    density = 10 # process every 10 lines
    resolution = 0.01 # mm/pixel
    dofilter =1 # filter the image
    notes = 8 # notes per octave
    maxscale = 8 #Max scale as inverse fraction of data length
    verbose = 1 # print stuff to screen
    dgs_stats = DGS.dgs(image_file, density, resolution, dofilter, maxscale, notes, verbose)

 REQUIRED INPUTS:
 simply a single file path
 
 OPTIONAL INPUTS [default values][range of acceptable values]
 density = process every *density* lines of image [10][1 - 100]
 resolution = spatial resolution of image in mm/pixel [1][>0]
 dofilter = spatial resolution of image in mm/pixel [1][0 or 1]
 notes = notes per octave to consider in continuous wavelet transform [8][1 - 8]
 maxscale = maximum scale (pixels) as an inverse function of data (image row) length [8][2 - 40]
 verbose = if 1, print stuff to screen [0][0 or 1]

OUTPUT:
A dictionary objects containing the following key/value pairs:
* mean grain size: arithmetic mean grain size
* grain size sorting: arithmetic standard deviation of grain sizes
* grain size skewness: arithmetic skewness of grain size-distribution
* grain size kurtosis: arithmetic kurtosis of grain-size distribution
* percentiles: 5th, 10th, 16th, 25th, 50th, 75th, 84th, 90th, and 95th percentile of the cumulative grain size (% less than) particle size distribution
* grain size frequencies: the normalised frequencies associated with 'grain size bins'
* grain size bins: grain size values at which the distribution is evaluated


PROCESSING NOTES:
Note that the larger the density parameter, the longer the execution time. 

:license:
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

import numpy as np
import sys #, getopt, os, glob
#from PIL.Image import open as imopen
import cwt
import sgolay
#from scipy.misc import imread as imopen
from imread import imread

from scipy.ndimage.interpolation import zoom

# suppress divide and invalid warnings
np.seterr(divide='ignore')
np.seterr(invalid='ignore')

import warnings
warnings.filterwarnings("ignore")

# =========================================================
def iseven(n):
   """Return true if n is even."""
   return n%2==0

# =========================================================
def isodd(n):
   """Return true if n is odd."""   
   return not iseven(n)

# =========================================================
def rescale(dat,mn,mx):
    """
    rescales an input dat between mn and mx
    """
    m = min(dat.flatten())
    M = max(dat.flatten())
    return (mx-mn)*(dat-m)/(M-m)+mn

# =========================================================
def get_me(useregion, maxscale, notes, density): #, mult):
   complete=0
   while complete==0:
      try:
         dat = cwt.Cwt(np.asarray(useregion,'int8'), maxscale, notes, density) #, mult)
         if 'dat' in locals(): 
            complete=1
      except:
         density = density +1

   return dat.getvar(), (np.pi/2)*dat.getscales()

# =========================================================
def filter_me(region):

   region = zoom(region, 0.5)

   nx, ny = np.shape(region)
   mn = min(nx,ny)

   if isodd(mn/4):
        window_size = (int(mn/4))
   else:
        window_size = (int(mn/4))-1

   if iseven(window_size):
      window_size = window_size+1

   Zf = sgolay.sgolay2d( region, window_size, order=3).getdata()
   # rescale filtered image to full 8-bit range
   useregion = rescale(zoom(region-Zf[:nx,:ny],2),0,255)

   return useregion

# =========================================================
# =========================================================
def dgs(image, density=10, resolution=1, dofilter=1, maxscale=8, notes=8, verbose=0):

   if verbose==1:
      print "==========================================="
      print "======DIGITAL GRAIN SIZE: WAVELET=========="
      print "==========================================="
      print "=CALCULATE GRAIN SIZE-DISTRIBUTION FROM AN="
      print "====IMAGE OF SEDIMENT/GRANULAR MATERIAL===="
      print "==========================================="
      print "======A PROGRAM BY DANIEL BUSCOMBE========="
      print "========USGS, FLAGSTAFF, ARIZONA==========="
      print "========REVISION 3.0.3, MAR 2016==========="
      print "==========================================="

   # exit program if no input folder given
   if not image:
      print 'An image file is required!!!!!!'
      sys.exit(2)

   # print given arguments to screen and convert data type where necessary
   if image:
      print 'Input image is ', image

   if density:
      density = np.asarray(density,int)
      print 'Every %s rows will be processed' % (str(density))

   if resolution:
      resolution = np.asarray(resolution,float)
      print 'Resolution is ', str(resolution)

   if dofilter:
      dofilter = np.asarray(dofilter,int)
      if dofilter==1:
         print 'Image will be filtered'
      else:
         print 'Image will not be filtered'

   if maxscale:
      maxscale = np.asarray(maxscale,int)
      print 'Max scale as inverse fraction of data length: %s' % (str(maxscale))

   if notes:
      notes = np.asarray(notes,int)
      print 'Analysis of %s sub-octaves per octave' % (str(notes))


   # ======= stage 1 ==========================
   # read image
   if verbose==1:
      print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
      print "Processing image %s" % (image)   
   try:
       #im = imopen(image, flatten=1).astype('uint8')#.convert("L")
       #im = imread.imload(image, as_grey=True).astype('uint8')
       im = imread(image)[:,:,:3] # read image up to 3 layers
       im = np.squeeze(im) # squeeze singleton dimensions
       if len(np.shape(im))==3: # if rgb, convert to grey
          im = (0.299 * im[:,:,0] + 0.5870*im[:,:,1] + 0.114*im[:,:,2]).astype('uint8')

       nx,ny = np.shape(im)
       if nx>ny:
          im=im.T

   except: # IOError:
       print 'cannot open', image
       sys.exit(2)
       #im = imread(image)

       #nx,ny = np.shape(im)
       #if nx>ny:
       #   im=im.T

   # convert to numpy array
   region = np.array(im)
   #nx, ny = np.shape(region)
   #mn = min(nx,ny)

   # ======= stage 2 ==========================
   # if requested, call sgolay to filter image
   if dofilter==1:
      useregion = filter_me(region) #, mn, nx, ny)

   else: #no filtering
      useregion = rescale(region,0,255)

   # ======= stage 3 ==========================
   # call cwt to get particle size distribution

   while (np.shape(useregion)[0] / density) > 100:
      density = density+1

   d, scales = get_me(useregion, maxscale, notes, density) #mult

   d = d/np.sum(d)
   d = d/(scales**0.5)
   d = d/np.sum(d)

   # ======= stage 4 ==========================
   # trim particle size bins
   index = np.nonzero(scales<ny/4)
   scales = scales[index]
   d = d[index]
   d = d/np.sum(d)

   index = np.nonzero(scales>np.pi*2)
   scales = scales[index]
   d = d[index]
   d = d/np.sum(d)
     
   n = np.r_[0:len(scales)]-(len(scales)-1)/2
   d = d*np.exp(-(0.5)*((np.pi/2)*n/((len(scales)-1)/2))**2)
   d = d/np.sum(d)   

   # get real scales by multiplying by resolution (mm/pixel)
   scales = scales*resolution

   # area-by-number to volume-by-number
   x = -1 # conversion constant
   r_v = (d*scales**x) / np.sum(d*scales**x) #volume-by-weight proportion

   # ======= stage 5 ==========================
   # calc particle size stats
   mnsz = np.sum(d*scales)
   if verbose==1:
      print "mean size = ", mnsz 

   srt = np.sqrt(np.sum(d*((scales-mnsz)**2)))
   if verbose==1:
      print "stdev = ",srt 

   sk = (sum(d*((scales-mnsz)**3)))/(100*srt**3)
   if verbose==1:
      print "skewness = ",sk

   kurt = (sum(d*((scales-mnsz)**4)))/(100*srt**4)
   if verbose==1:
      print "kurtosis = ",kurt

   pd = np.interp([.05,.1,.16,.25,.5,.75,.84,.9,.95],np.hstack((0,np.cumsum(d))), np.hstack((0,scales)) )

   # ======= stage 6 ==========================
   # return a dict object of stats
   return {'mean grain size': mnsz, 'grain size sorting': srt, 'grain size skewness': sk, 'grain size kurtosis': kurt, 'percentiles': [.05,.1,.16,.25,.5,.75,.84,.9,.95], 'percentile_values': pd, 'grain size frequencies': d, 'grain size bins': scales}


# =========================================================
# =========================================================
if __name__ == '__main__':

   dgs(image, density=10, resolution=1, dofilter=1, maxscale=8, notes=8, verbose=0)


