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

from dgs import *
import os, glob
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('fivethirtyeight')

#========================================
# single image
def dotest1(image, with_plot=False):
   # if this is 1, it means "give me the results in pixels - I'll apply my own scaling"
   # otherwise, it is mm/pixel (if you want your results in mm) or um/pixel for microns
   resolution = 1 #.04

   #the maximum scale (grain size) considered by the wavelet is the horizontal width dimension divided by this number
   #so if your image is 1000 pixels wide and maxscale=4, only grains up to 1000/4 = 250 pixels are considered
   maxscale=10 

   # if 1, prints grain size statistics to screen
   verbose=1

   #this is the area to volume conversion coefficient. See Cuttler et al (provided)
   #you could also use it as an empirical tuning coefficient against field data (recommended)
   x = 0

   # 0 means do not apply denoising filter
   filter = 0

   #I recommend you compute in pixels (resolution=1) then apply your resolution scaling afterwards
   data_out = dgs(image, resolution, maxscale, verbose, x, filter)

   ## parse out dict into three separate dictionaries
   stats = dict(list(data_out.items())[:4])
   percentiles = dict(list(data_out.items())[4:6])
   freqs_bins = dict(list(data_out.items())[6:])

   if resolution!=1:
       freqs_bins['grain size bins']*=resolution
       percentiles['percentile_values']*=resolution

       for k in stats.keys():
           stats[k] = stats[k]*resolution

   # write each to csv file
   pd.DataFrame.from_dict(stats.items()).to_csv('demo_results/'+image.split(os.sep)[-1]+'_stats.csv')
   pd.DataFrame.from_dict(percentiles).to_csv('demo_results/'+image.split(os.sep)[-1]+'_percentiles.csv')
   pd.DataFrame.from_dict(freqs_bins).to_csv('demo_results/'+image.split(os.sep)[-1]+'_freqs_bins.csv')

   if with_plot == True:
       #do stufff
       plt.plot(freqs_bins['grain size bins'], freqs_bins['grain size frequencies'],'k', lw=2, label=image)
       plt.legend()
       if resolution!=1:
           plt.xlabel('Grain Size (mm or units provided)')
       else:
           plt.xlabel('Grain Size (pixels)')
       plt.ylabel('Frequency')
       #plt.show()
       plt.savefig('demo_results/1image_psd.png', dpi=300, bbox_inches='tight')
       plt.close('all')


#========================================
## folder of images
def dotest_batch(folder, set=1, with_plot=False):

   if set==1:
      files = glob.glob(folder+os.sep+'IMG*.JPG')
      files = [f for f in files if f.endswith('.JPG')]
      resolution = 1#.04
      maxscale=10 
   elif set==2:
      files = glob.glob(folder+os.sep+'IMG*.jpg')
      files = [f for f in files if f.endswith('.jpg')]
      resolution = 1#.04
      maxscale=10 
   elif set==3:
      files = glob.glob(folder+os.sep+'*.tif')
      resolution =1# .04
      maxscale=10

   # if this is 1, it means "give me the results in pixels - I'll apply my own scaling"
   # otherwise, it is mm/pixel (if you want your results in mm) or um/pixel for microns

   #the maximum scale (grain size) considered by the wavelet is the horizontal width dimension divided by this number
   #so if your image is 1000 pixels wide and maxscale=4, only grains up to 1000/4 = 250 pixels are considered

   # if 1, prints grain size statistics to screen
   verbose=0

   #this is the area to volume conversion coefficient. See Cuttler et al (provided)
   #you could also use it as an empirical tuning coefficient against field data (recommended)
   x = 0

   # 1 means apply denoising filter
   filter = 1

   ALL_RES = []
   for f in tqdm(files): #tqdm gives you a progress bar
      data_out = dgs(f, resolution, maxscale, verbose, x, filter)
      ALL_RES.append(data_out)


   ## parse out dict into three separate dictionaries
   S = {}; P = {}; F = {}
   counter = 0
   for data_out in ALL_RES:
      stats = dict(list(data_out.items())[:4])
      percentiles = dict(list(data_out.items())[4:6])
      freqs_bins = dict(list(data_out.items())[6:])

      if resolution!=1:
         freqs_bins['grain size bins']*=resolution
         percentiles['percentile_values']*=resolution

         for k in stats.keys():
             stats[k] = stats[k]*resolution

      S[files[counter]] = stats.items()
      P[files[counter]] = percentiles
      F[files[counter]] = freqs_bins
      counter += 1

   # convert into stats (rows) versus images (columns)
   tmp = list(S.keys())
   d = {tmp[0]: [k[1] for k in list(S[tmp[0]])]}
   for k in range(1,len(tmp)):
       d.update( {tmp[k]: [k[1] for k in list(S[tmp[k]])]} )

   pd.DataFrame(data=d, index = ['mean grain size', 'grain size sorting', 'grain size skewness', 'grain size kurtosis']).to_csv('demo_results/stats_batch.csv')

   # convert into percentiles (rows) versus images (columns)
   tmp = list(P.keys())
   d = {tmp[0]: P[tmp[0]]['percentile_values']}
   for k in range(1,len(tmp)):
       d.update( {tmp[k]: P[tmp[k]]['percentile_values'] } )

   pd.DataFrame(data=d, index = P[tmp[0]]['percentiles']).to_csv('demo_results/percentiles_batch.csv')

   # write each to csv file
   # pd.DataFrame.from_dict(S).to_csv('demo_results/stats_batch.csv')
   # pd.DataFrame.from_dict(P).to_csv('demo_results/percentiles_batch.csv')
   pd.DataFrame.from_dict(F).to_csv('demo_results/freqs_bins_batch.csv')

   if with_plot == True:
       counter = 0
       cols = ['r','g','b','m','c','k','y'][:len(F)]
       for f in F:
          plt.plot(F[f]['grain size bins'], F[f]['grain size frequencies'],cols[counter], lw=2, label=files[counter])
          counter += 1
       plt.legend()

       if resolution!=1:
           plt.xlabel('Grain Size (mm or units provided)')
       else:
           plt.xlabel('Grain Size (pixels)')

       #plt.xlabel('Grain Size (pixels)')
       plt.ylabel('Frequency')
       #plt.show()
       if set==1:
          plt.savefig('demo_results/batch_sand_3images_psd.png', dpi=300, bbox_inches='tight')
       elif set==2:
          plt.savefig('demo_results/batch_sand_6images_psd.png', dpi=300, bbox_inches='tight')
       elif set==3:
          plt.savefig('demo_results/batch_gravel_6images_psd.png', dpi=300, bbox_inches='tight')
       plt.close('all')


#====================================
if __name__ == '__main__':

   #image= 'data'+os.sep+'IMG_0249.JPG'  #finr
   image= 'data'+os.sep+'IMG_0229.JPG'  #medium
   #image= 'data'+os.sep+'IMG_0202.JPG'   #big

   # one image, no plot
   dotest1(image)

   image= 'data'+os.sep+'IMG_0202.JPG'   #big

   # new image, with plot
   dotest1(image, True)

   # all images in data folder, with plot
   folder = 'data'
   set = 1
   dotest_batch(folder, set, True)

   set = 2
   dotest_batch(folder, set, True)

   set = 3
   dotest_batch(folder, set, True)
