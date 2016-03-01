### About

![alt tag](http://dbuscombe-usgs.github.io/figs/2013-02-24-dgs/nicegrains1.jpg)

pyDGS - a Python framework for wavelet-based digital grain size analysis

pyDGS is an open-source project dedicated to provide a Python framework to compute estimates of grain size distribution  using the continuous wavelet transform method of Buscombe (2013) from an image of sediment where grains are clearly resolved. DOES NOT REQUIRE CALIBRATION

This program implements the algorithm of:

Buscombe, D. (2013) Transferable Wavelet Method for Grain-Size Distribution from Images of Sediment Surfaces and Thin Sections, and Other Natural Granular Patterns. Sedimentology 60, 1709-1732

http://dbuscombe-usgs.github.io/docs/Buscombe2013_Sedimentology_sed12049.pdf

### install:
    python setup.py install
    sudo python setup.py install
    pip install pyDGS
    
### test:
    python -c "import DGS; DGS.test.dotest()"

 REQUIRED INPUTS:
 image name e.g. '/home/sed_images/my_image.png'
 
 OPTIONAL INPUTS [default values][range of acceptable values]
 density = process every density lines of image [10][1 - 100]
 resolution = spatial resolution of image in mm/pixel [1][>0]
 dofilter = spatial resolution of image in mm/pixel [1][0 or 1]
 notes = notes per octave to consider in continuous wavelet transform [8][1 - 8]
 maxscale = maximum scale (pixels) as an inverse function of data (image row) length [8][2 - 40]
 doplot = 0=no, 1=yes [0][0 or 1]

OUTPUT FOR A SINGLE IMAGE FILE:
A dictionary objects containing the following key/value pairs:
* mean grain size: arithmetic mean grain size
* grain size sorting: arithmetic standard deviation of grain sizes
* grain size skewness: arithmetic skewness of grain size-distribution
* grain size kurtosis: arithmetic kurtosis of grain-size distribution
* percentiles: 5th, 10th, 16th, 25th, 50th, 75th, 84th, 90th, and 95th percentile of the cumulative grain size (% less than) particle size distribution
* grain size frequencies: the normalised frequencies associated with 'grain size bins'
* grain size bins: grain size values at which the distribution is evaluated


### processing example on 1 image:
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

### license:
    GNU Lesser General Public License, Version 3
    (http://www.gnu.org/copyleft/lesser.html)
    
    This software is in the public domain because it contains materials that
    originally came from the United States Geological Survey, an agency of the
    United States Department of Interior. For more information, 
    see the official USGS copyright policy at
    http://www.usgs.gov/visual-id/credit_usgs.html#copyright
    Any use of trade, product, or firm names is for descriptive purposes only 
    and does not imply endorsement by the U.S. government.
    
### Note for Windows Users

I recommend the Anaconda python distribution for Windows which includes all of the library dependencies required to run this program. Anaconda comes with a variety of IDEs and is pretty easy to use. To run the test images, launch the Anaconda command terminal and type:

```
pip install pyDGS
python -c "import DGS; DGS.test.dotest()"
```

### Contributing & Credits

This program implements the algorithm of 
Buscombe, D. (2013) Transferable Wavelet Method for Grain-Size Distribution from Images of Sediment Surfaces and Thin Sections, and Other Natural Granular Patterns, Sedimentology 60, 1709 - 1732

 Author:  Daniel Buscombe
           Grand Canyon Monitoring and Research Center
           United States Geological Survey
           Flagstaff, AZ 86001
           dbuscombe@usgs.gov
 Revision Mar 1, 2015
 First Revision January 18 2013

For more information visit https://github.com/dbuscombe-usgs/pyDGS

Please contact:
dbuscombe@usgs.gov

to report bugs and discuss the code, algorithm, collaborations

For the latest code version please visit:
https://github.com/dbuscombe-usgs

See also the project blog: 
http://dbuscombe-usgs.github.com/

Please download, try, report bugs, fork, modify, evaluate, discuss. Thanks for stopping by!
