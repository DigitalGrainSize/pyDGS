
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
   resolution = 1

   #the maximum scale (grain size) considered by the wavelet is the horizontal width dimension divided by this number
   #so if your image is 1000 pixels wide and maxscale=4, only grains up to 1000/4 = 250 pixels are considered
   maxscale=4

   # if 1, prints grain size statistics to screen
   verbose=1

   #this is the area to volume conversion coefficient. See Cuttler et al (provided)
   #you could also use it as an empirical tuning coefficient against field data (recommended)
   x = 0

   data_out = dgs(image, resolution, maxscale, verbose, x)

   ## parse out dict into three separate dictionaries
   stats = dict(list(data_out.items())[:4])
   percentiles = dict(list(data_out.items())[4:6])
   freqs_bins = dict(list(data_out.items())[6:])

   # write each to csv file
   pd.DataFrame.from_dict(stats.items()).to_csv('demo_results/stats.csv')
   pd.DataFrame.from_dict(percentiles).to_csv('demo_results/percentiles.csv')
   pd.DataFrame.from_dict(freqs_bins).to_csv('demo_results/freqs_bins.csv')

   if with_plot == True:
       #do stufff
       plt.plot(freqs_bins['grain size bins'], freqs_bins['grain size frequencies'],'k', lw=2, label=image)
       plt.legend()
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
   elif set==2:
      files = glob.glob(folder+os.sep+'IMG*.jpg')
      files = [f for f in files if f.endswith('.jpg')]
   elif set==3:
      files = glob.glob(folder+os.sep+'*.tif')

   # if this is 1, it means "give me the results in pixels - I'll apply my own scaling"
   # otherwise, it is mm/pixel (if you want your results in mm) or um/pixel for microns
   resolution = 1

   #the maximum scale (grain size) considered by the wavelet is the horizontal width dimension divided by this number
   #so if your image is 1000 pixels wide and maxscale=4, only grains up to 1000/4 = 250 pixels are considered
   maxscale=4

   # if 1, prints grain size statistics to screen
   verbose=0

   #this is the area to volume conversion coefficient. See Cuttler et al (provided)
   #you could also use it as an empirical tuning coefficient against field data (recommended)
   x = 0

   ALL_RES = []
   for f in tqdm(files): #tqdm gives you a progress bar
      data_out = dgs(f, resolution, maxscale, verbose, x)
      ALL_RES.append(data_out)


   ## parse out dict into three separate dictionaries
   S = {}; P = {}; F = {}
   counter = 0
   for data_out in ALL_RES:
      stats = dict(list(data_out.items())[:4])
      percentiles = dict(list(data_out.items())[4:6])
      freqs_bins = dict(list(data_out.items())[6:])
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
       plt.xlabel('Grain Size (pixels)')
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
