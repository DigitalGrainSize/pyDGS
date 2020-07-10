
## About

pydgs is an open-source project dedicated to provide a Python framework to compute estimates of grain size distribution  using the continuous wavelet transform method of Buscombe (2013) from an image of sediment where grains are clearly resolved. It doesn't require calibration, but does require scaling. It works best for well-sorted sands and gravels. It doesn't work so well for mixtures or subpixel grains - you might have better luck with my other program [SediNet](https://github.com/MARDAScience/SediNet)

This program implements the algorithm of:

Buscombe, D. (2013) Transferable Wavelet Method for Grain-Size Distribution from Images of Sediment Surfaces and Thin Sections, and Other Natural Granular Patterns. Sedimentology 60, 1709-1732]. See also Cuttler et al., (2017) for details on the implementation of the area-by-number to volume-by-number conversion. Pdfs are provided in `docs`.


## Install:

download the code and change directory

```
git clone --depth 1 https://github.com/dbuscombe-usgs/pyDGS.git
cd pyDGS
```

install the provided conda environment

```
conda env create -f conda_env/pydgs.yml
conda activate pydgs
```

Run tests:
```
python test.py
```

Test results are written to `demo_results`. Imagery is in `data`

Adapt `test.py` to your own needs, to analyze your own imagery


### REQUIRED INPUTS:

 image name e.g. `'/home/sed_images/my_image.png'`

### OPTIONAL INPUTS

 * `resolution` = spatial resolution of image in mm/pixel `[1][>0]`. For results in pixels, use resolution = 1
 * `verbose` = if 1, prints grain size statistics to screen `[0][0 or 1]`
 * `x` = area-by-number to volume-by-number conversion `[0] [-1 - +1]`
 * `maxscale` =  the maximum scale (grain size) considered by the wavelet is the horizontal width dimension divided by this number. So if your image is 1000 pixels wide and maxscale=4, only grains up to 1000/4 = 250 pixels are considered

See also Cuttler et al., 2017 (in `docs`) for details on the implementation of the area-by-number to volume-by-number conversion. You could also use it as an empirical tuning coefficient against field data (recommended)


### OUTPUT FOR A SINGLE IMAGE FILE:

A dictionary objects containing the following key/value pairs:
* mean grain size: arithmetic mean grain size
* grain size sorting: arithmetic standard deviation of grain sizes
* grain size skewness: arithmetic skewness of grain size-distribution
* grain size kurtosis: arithmetic kurtosis of grain-size distribution
* percentiles: 5th, 10th, 16th, 25th, 50th, 75th, 84th, 90th, and 95th percentile of the cumulative grain size (% less than) particle size distribution
* grain size frequencies: the normalised frequencies associated with 'grain size bins'
* grain size bins: grain size values at which the distribution is evaluated


## Contributing & Credits

This program implements the algorithm of
Buscombe, D. (2013) Transferable Wavelet Method for Grain-Size Distribution from Images of Sediment Surfaces and Thin Sections, and Other Natural Granular Patterns, Sedimentology 60, 1709 - 1732

> **Author**:  Daniel Buscombe  
>          Marda Science, LLC
>          (formerly of NAU and USGS and Plymouth University)
>          Flagstaff, AZ   

 Revision: July 2020
 First Revision: January 18 2013

For more information visit <https://github.com/dbuscombe-usgs/pyDGS>


## Major changes in version 4 (July 2020)

* switched to using the `pywavelets` package ([see here](https://pywavelets.readthedocs.io/en/latest/)) for wavelet calculations, which removes many issues with my previous cython translation, including no need to specify parameters
* more test data
* improved image filtering, which now occurs as standard
* improved test script to provide example usage
* simpler to install. no code compilation, no pip support, no `setup.py` script
