
## About

pyDGS is an open-source project dedicated to provide a Python framework to compute estimates of grain size distribution  using the continuous wavelet transform method of Buscombe (2013) from an image of sediment where grains are clearly resolved. DOES NOT REQUIRE CALIBRATION

This program implements the algorithm of:

[Buscombe, D. (2013) Transferable Wavelet Method for Grain-Size Distribution from Images of Sediment Surfaces and Thin Sections, and Other Natural Granular Patterns. Sedimentology 60, 1709-1732](https://www.danielbuscombe.com/s/Buscombe_2013_sedimentology_101111-sed12049.pdf)

See also [Cuttler et al., 2017](https://www.danielbuscombe.com/s/Cuttler_et_al-2017-Sedimentology.pdf) for details on the implementation of the area-by-number to volume-by-number conversion


## Install:

``` bash
pip install git+https://github.com/dbuscombe-usgs/pyDGS.git
```

Or, install the provided conda environment

``` bash
conda env create -f binder\environment.yml
conda activate DGS
cd ..
pip install git+https://github.com/dbuscombe-usgs/pyDGS.git
```

If you intend to run the jupyter examples (see below), you may also need to install a jupyter kernel associated with the conda environment

```
python -m ipykernel install --user --name DGS --display-name "Python (dgs)"
```

## Run the jupyter examples

```jupyter notebook```

then head to the browser, navigate to notebooks folder, and launch a notebook

## Run on binder

https://mybinder.org/v2/gh/dbuscombe-usgs/pyDGS/master


## Run on Google Cloud Platform

First, follow instructions here for how to set up an instance to run in GCP. Make sure to set a static IP address, as per the instructions, and make a note of that because you'll need it later

Then open a shell into the VM and set it up to

```
ssh-keygen -t rsa -b 4096 -C "yourname@youremail.com"

eval "$(ssh-agent -s)"

ssh-add ~/.ssh/id_rsa

cat ~/.ssh/id_rsa.pub
```

Then copy the key into your github profile keys. For more information about how to do that, see here. xclip likely won't work, but you can simply copy (Ctrl-C) the text printed to screen

You will be cloning your fork of the main repo, so replace YOURUSERNAME in the below code to clone the repo and set up a conda environment to run in

``` bash
cd ..
pip install git+https://github.com/YOURUSERNAME/pyDGS.git
```

Now you can run pyDGS on the cloud.

To run the jupyter notebooks, run the following command to run the jupyter notebook server

```
cd pyDGS/notebooks
jupyter notebook --NotebookApp.iopub_data_rate_limit=10000000
```

The jupyterlab server will be displayed at

```
http://IP:8888
```

where IP is the static IP of the VM that you noted earlier.



## Test:

You can run pyDGS tests using the following script:

``` bash
python -c "import DGS; DGS.test.dotest()"
```

### REQUIRED INPUTS:

 image name e.g. `'/home/sed_images/my_image.png'`

### OPTIONAL INPUTS 

#### ***`[default values][range of acceptable values]`***

 * **density** = process every *density* lines of image `[10][1 - 100]`
 * **resolution** = spatial resolution of image in mm/pixel `[1][>0]`
 * **dofilter** = spatial resolution of image in mm/pixel `[1][0 or 1]`
 * **notes** = notes per octave to consider in continuous wavelet transform `[8][1 - 8]`
 * **maxscale** = maximum scale (pixels) as an inverse function of data (image row) length `[8][2 - 40]`
 * **verbose** = if 1, print stuff to screen `[0][0 or 1]`
 * **x** = area-by-number to volume-by-number conversion `[0] [-1 - +1]`

See also [Cuttler et al., 2017](https://www.danielbuscombe.com/s/Cuttler_et_al-2017-Sedimentology.pdf) for details on the implementation of the area-by-number to volume-by-number conversion


### OUTPUT FOR A SINGLE IMAGE FILE:

A dictionary objects containing the following key/value pairs:
* mean grain size: arithmetic mean grain size
* grain size sorting: arithmetic standard deviation of grain sizes
* grain size skewness: arithmetic skewness of grain size-distribution
* grain size kurtosis: arithmetic kurtosis of grain-size distribution
* percentiles: 5th, 10th, 16th, 25th, 50th, 75th, 84th, 90th, and 95th percentile of the cumulative grain size (% less than) particle size distribution
* grain size frequencies: the normalised frequencies associated with 'grain size bins'
* grain size bins: grain size values at which the distribution is evaluated


### Processing example on 1 image:

From a terminal with pyDGS already installed using `pip`:

``` bash
$ python

Python 3.6.5 (default, Apr 11 2018, 10:42:01)
[GCC 4.2.1 Compatible Apple LLVM 9.1.0 (clang-902.0.39.1)] on darwin
Type "help", "copyright", "credits" or "license" for more information.
>>>
```

Then, once you're in the python terminal (or as a separate .py file):

``` python
import DGS

image_file = '/home/sed_images/my_image.png'

density = 10 # process every 10 lines
resolution = 0.01 # mm/pixel
dofilter =1 # filter the image
notes = 16 # notes per octave
maxscale = 2 #Max scale as inverse fraction of data length
verbose = 1 # print stuff to screen
x = 1
dgs_stats = DGS.dgs(image_file, density, resolution, dofilter, maxscale, notes, verbose, x)
```

#### REQUIRED INPUTS:

 simply a single file path

### OPTIONAL INPUTS 

#### ***[default values][range of acceptable values]***

* **density** = process every *density* lines of image [10][1 - 100]
* **resolution** = spatial resolution of image in mm/pixel [1][>0]
* **dofilter** = spatial resolution of image in mm/pixel [1][0 or 1]
* **notes** = notes per octave to consider in continuous wavelet transform [16][1 - 8]
* **maxscale** = maximum scale (pixels) as an inverse function of data (image row) length [2][2 - 40]
* **verbose** = if 1, print stuff to screen [0][0 or 1]
* **x** = calibration coefficient (nominally, area-by-number to volume-by-number conversion) [1] [-3 - 3]

### OUTPUT:

A dictionary objects containing the following key/value pairs:

* **mean grain size**: arithmetic mean grain size
* **grain size sorting**: arithmetic standard deviation of grain sizes
* **grain size skewness**: arithmetic skewness of grain size-distribution
* **grain size kurtosis**: arithmetic kurtosis of grain-size distribution
* **percentiles**: 5th, 10th, 16th, 25th, 50th, 75th, 84th, 90th, and 95th percentile of the cumulative grain size (% less than) particle size distribution
* **grain size frequencies**: the normalised frequencies associated with 'grain size bins'
* **grain size bins**: grain size values at which the distribution is evaluated

### PROCESSING NOTES:

Note that the larger the density parameter, the longer the execution time. 

## License:
    GNU Lesser General Public License, Version 3
    (http://www.gnu.org/copyleft/lesser.html)
    
    This software is in the public domain because it contains materials that
    originally came from the United States Geological Survey, an agency of the
    United States Department of Interior. For more information, 
    see the official USGS copyright policy at
    http://www.usgs.gov/visual-id/credit_usgs.html#copyright
    Any use of trade, product, or firm names is for descriptive purposes only 
    and does not imply endorsement by the U.S. government.


## Contributing & Credits

This program implements the algorithm of 
Buscombe, D. (2013) Transferable Wavelet Method for Grain-Size Distribution from Images of Sediment Surfaces and Thin Sections, and Other Natural Granular Patterns, Sedimentology 60, 1709 - 1732

> **Author**:  Daniel Buscombe  
>          Northern Arizona University  
>          Flagstaff, AZ 86001  
>          daniel.buscombe@nau.edu

 Revision: July 2019
 First Revision: January 18 2013

For more information visit <https://github.com/dbuscombe-usgs/pyDGS>

<https://www.danielbuscombe.com/s/Buscombe_2013_sedimentology_101111-sed12049.pdf>

Please contact:
<daniel.buscombe@nau.edu>

to report bugs and discuss the code, algorithm, collaborations

> DO NOT EMAIL dbuscombe@usgs.gov - that email is never checked!!!

For the latest code version please visit:
<https://github.com/dbuscombe-usgs>

See also the project blog: 
<http://dbuscombe-usgs.github.com/>

Please download, try, report bugs, fork, modify, evaluate, discuss. Thanks for stopping by!
