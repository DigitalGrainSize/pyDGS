# Written by Dr Daniel Buscombe, Marda Science LLC
#
# MIT License
#
# Copyright (c) 2020, Marda Science LLC
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
import sys, getopt
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('fivethirtyeight')
from tkinter import Tk
from tkinter.filedialog import askopenfilename, askdirectory
from datetime import datetime

#================================================================
def do_dgs(resolution, maxscale, verbose, files):

   ALL_RES = []
   for f in tqdm(files): #tqdm gives you a progress bar
      data_out = dgs(f, 1, maxscale, verbose, x)
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

   timestr = datetime.now().strftime("%Y-%m-%d-%H-%M")

   pd.DataFrame(data=d, index = ['mean grain size', 'grain size sorting', 'grain size skewness', 'grain size kurtosis']).to_csv('demo_results/stats_batch_'+timestr+'.csv')

   # convert into percentiles (rows) versus images (columns)
   tmp = list(P.keys())
   d = {tmp[0]: P[tmp[0]]['percentile_values']}
   for k in range(1,len(tmp)):
       d.update( {tmp[k]: P[tmp[k]]['percentile_values'] } )

   pd.DataFrame(data=d, index = P[tmp[0]]['percentiles']).to_csv('demo_results/percentiles_batch_'+timestr+'.csv')

   # write each to csv file
   # pd.DataFrame.from_dict(S).to_csv('demo_results/stats_batch.csv')
   # pd.DataFrame.from_dict(P).to_csv('demo_results/percentiles_batch.csv')
   pd.DataFrame.from_dict(F).to_csv('demo_results/freqs_bins_batch_'+timestr+'.csv')

   counter = 0
   cols = ['r','g','b','m','c','k','y'][:len(F)]
   for f in F:
      plt.plot(F[f]['grain size bins'], F[f]['grain size frequencies'],cols[counter], lw=2, label=files[counter].split(os.sep)[-1])
      counter += 1
   plt.legend(fontsize=6)

   if resolution!=1:
       plt.xlabel('Grain Size (mm)')
   else:
       plt.xlabel('Grain Size (pixels)')

   #plt.xlabel('Grain Size (pixels)')
   plt.ylabel('Frequency')
   #plt.show()
   plt.savefig('demo_results/batch_psd_'+timestr+'.png', dpi=300, bbox_inches='tight')
   plt.close('all')

#====================================
if __name__ == '__main__':

    argv = sys.argv[1:]
    try:
        opts, args = getopt.getopt(argv,"h:r:m:x:")
    except getopt.GetoptError:
        print('======================================')
        print('python run_dgs.py') #
        print('python run_dgs.py {-r resolution in mm per pixel (float)} {-m maxscale *see below (integer)} {-x "x" parameter **see below (float) }') #
        print('*the maximum scale (grain size) considered by the wavelet is the horizontal width dimension divided by this number') #
        print('so if your image is 2000 pixels wide and maxscale=8, only grains up to 2000/8 = 250 pixels are considered')
        print('**this is the area to volume conversion coefficient. See Cuttler et al (provided)')
        print('you could also use it as an empirical tuning coefficient against field data (recommended)')
        print('======================================')

        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('======================================')
            print('python run_dgs.py') #
            print('python run_dgs.py {-r resolution in mm per pixel (float)} {-m maxscale *see below (integer)} {-x "x" parameter **see below (float) }') #
            print('*the maximum scale (grain size) considered by the wavelet is the horizontal width dimension divided by this number') #
            print('so if your image is 2000 pixels wide and maxscale=8, only grains up to 2000/8 = 250 pixels are considered')
            print('**this is the area to volume conversion coefficient. See Cuttler et al (provided)')
            print('you could also use it as an empirical tuning coefficient against field data (recommended)')
            print('======================================')
            print('======================================')
            print('Example usage: python run_dgs.py -r 0.04')
            print('Example usage: python run_dgs.py -m 20')
            print('Example usage: python run_dgs.py -r 0.04 -m 10 -x 0.5')
            print('Example usage: python run_dgs.py -r 0.04 -m 20 -x -0.1')
            print('Example usage: python run_dgs.py -x -0.5')
            print('======================================')
            sys.exit()
        elif opt in ("-r"):
            resolution = arg
            resolution = float(resolution)
        elif opt in ("-r"):
            maxscale = arg
            maxscale = int(maxscale)
        elif opt in ("-x"):
            x = arg
            x = float(x)

    if 'resolution' not in locals():
        resolution = 1
        print('Warning: no resolution in mm/px specified, using %i by default' % (resolution))
    if 'maxscale' not in locals():
        maxscale = 20
        print('Warning: specify a maxscale for best results, using %i by default' % (maxscale))
    if 'x' not in locals():
        x = 0.0
        print('Warning: specify "x" for best results, using %f by default' % (x))

    Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
    files = askopenfilename(title='Select image files', multiple=True, filetypes=[("Pick files","*.*")])

    # use verbose=1 for more output from dgs
    verbose=0

    # exit program if no input folder given
    if not files:
       print('Image files are required! ... program exiting')
       sys.exit(2)

    if resolution:
       resolution = np.asarray(resolution,float)
       print('Resolution is '+str(resolution))

    if maxscale:
       maxscale = np.asarray(maxscale,int)
       print('Max scale as inverse fraction of data length: '+str(maxscale))

    if x:
       x = np.asarray(x, float)
       print('Area to volume conversion constant = '+str(x))

    do_dgs(resolution, maxscale, verbose, files)

##
