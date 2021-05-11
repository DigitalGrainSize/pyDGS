
## About

pydgs is an open-source project dedicated to provide a Python framework to compute estimates of grain size distribution  using the continuous wavelet transform method of Buscombe (2013) from an image of sediment where grains are clearly resolved. It doesn't require calibration, but does require scaling. It works best for well-sorted sands and gravels. It doesn't work so well for mixtures or subpixel grains - you might have better luck with my other program [SediNet](https://github.com/MARDAScience/SediNet)

This program implements the algorithm of:

Buscombe, D. (2013) Transferable Wavelet Method for Grain-Size Distribution from Images of Sediment Surfaces and Thin Sections, and Other Natural Granular Patterns. Sedimentology 60, 1709-1732]. See also Cuttler et al., (2017) for details on the implementation of the area-by-number to volume-by-number conversion. Pdfs are provided in `docs`.

## Contents
* [Installation](#install)
* [Installation](#RUN_DGS.PY)
* [Inputs](#inputs)
* [Outputs](#outputs)
* [Credits](#ack)
* [Changelog](#change)


## <a name="install"></a>Install and test:

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

### <a name="RUN_DGS.PY"></a>RUN_DGS.PY:

Adapt `test.py` to your own needs, to analyze your own imagery

Or use the provided script

`python run_dgs.py`

Full syntax:

`python run_dgs.py {-r resolution in mm per pixel (float)} {-m maxscale *see below (integer)} {-x "x" parameter **see below (float) }`

*the maximum scale (grain size) considered by the wavelet is the horizontal width dimension divided by this number
so if your image is 2000 pixels wide and maxscale=8, only grains up to 2000/8 = 250 pixels are considered')

**this is the area to volume conversion coefficient. See Cuttler et al (provided)')
you could also use it as an empirical tuning coefficient against field data (recommended)')

Example usage

```
python run_dgs.py -r 0.04
python run_dgs.py -m 20
python run_dgs.py -r 0.04 -m 10 -x 0.5
python run_dgs.py -r 0.04 -m 20 -x -0.1
python run_dgs.py -x -0.5
```

### <a name="inputs"></a>REQUIRED INPUTS:

 image name e.g. `'/home/sed_images/my_image.png'`

### SUGGESTED INPUTS

 * `resolution` = spatial resolution of image in mm/pixel `[1][>0]`. For results in pixels, use resolution = 1
 * `x` = area-by-number to volume-by-number conversion `[0] [-1 - +1]`
 * `maxscale` =  the maximum scale (grain size) considered by the wavelet is the horizontal width dimension divided by this number. So if your image is 1000 pixels wide and maxscale=4, only grains up to 1000/4 = 250 pixels are considered

See also Cuttler et al., 2017 (in `docs`) for details on the implementation of the area-by-number to volume-by-number conversion. You could also use it as an empirical tuning coefficient against field data (recommended)


### <a name="outputs"></a> OUTPUT FOR A SINGLE IMAGE FILE:

A dictionary objects containing the following key/value pairs:
* mean grain size: arithmetic mean grain size
* grain size sorting: arithmetic standard deviation of grain sizes
* grain size skewness: arithmetic skewness of grain size-distribution
* grain size kurtosis: arithmetic kurtosis of grain-size distribution
* percentiles: 5th, 10th, 16th, 25th, 30th, 50th, 75th, 84th, 90th, and 95th percentile of the cumulative grain size (% less than) particle size distribution
* grain size frequencies: the normalised frequencies associated with 'grain size bins'
* grain size bins: grain size values at which the distribution is evaluated


## <a name="ack"></a> Contributing & Credits

This program implements the algorithm of
Buscombe, D. (2013) Transferable Wavelet Method for Grain-Size Distribution from Images of Sediment Surfaces and Thin Sections, and Other Natural Granular Patterns, Sedimentology 60, 1709 - 1732

> **Author**:  Daniel Buscombe  
>          Marda Science, LLC
>          (formerly of NAU and USGS and Plymouth University)
>          Flagstaff, AZ   

 Revision: May 2021
 First Revision: January 18 2013

For more information visit <https://github.com/dbuscombe-usgs/pyDGS>

## <a name="change"></a> CHANGELOG
### Major changes in version 4 (July 2020)

* switched to using the `pywavelets` package ([see here](https://pywavelets.readthedocs.io/en/latest/)) for wavelet calculations, which removes many issues with my previous cython translation, including no need to specify parameters
* more test data
* improved image filtering, which now occurs as standard
* improved test script to provide example usage
* simpler to install. no code compilation, no pip support, no `setup.py` script

### Changes in version 4.1 (May 2020)
* added `run_dgs.py` script utility to help implement the dgs module
* added image standardization - each input is now scaled by its mean and standard deviation such that it has zero mean and unit variance
