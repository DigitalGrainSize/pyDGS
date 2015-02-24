#!/usr/bin/env python
# -*- coding: utf-8 -*-
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
 Revision Feb 23, 2015
 First Revision January 18 2013   

For more information visit https://github.com/dbuscombe-usgs/pyDGS

:install:
    python setup.py install
    sudo python setup.py install
    pip install pyDGS
    
:test:
    python -c "import DGS; DGS.test.dotest()"

:processing example:
    python
    import DGS
    density = 10 # process every 10 lines
    res = 0.01 # mm/pixel
    doplot = 0 # don't make plots
    image_folder = '/home/sed_images'
    DGS.dgs(image_folder,density,doplot,res)
    image_file = '/home/sed_images/my_image.png'
    mnsz, srt, sk, kurt, pd = DGS.dgs(image_file,density,doplot,res)

 REQUIRED INPUTS:
 folder e.g. '/home/my_sediment_images'
 if 'pwd', then the present directory is analysed
 or simply a single file
 
 OPTIONAL INPUTS [default values]
 density = process every density lines of image [10]
 doplot = 0=no, 1=yes [0]
 resolution = spatial resolution of image in mm/pixel [1]

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
#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

import os, sys, glob
import inspect
from distutils.core import setup
from distutils.extension import Extension
#from setuptools import setup, Extension
import numpy as np

# Directory of the current file 
SETUP_DIRECTORY = os.path.dirname(os.path.abspath(inspect.getfile(
    inspect.currentframe())))

# Set this to True to enable building extensions using Cython.
# Set it to False to build extensions from the C file (that
# was previously created using Cython).
# Set it to 'auto' to build with Cython if available, otherwise
# from the C file.
USE_CYTHON = True

if USE_CYTHON:
   try:
      from Cython.Distutils import build_ext
   except:
      USE_CYTHON = False

# Read version from distmesh/__init__.py
with open(os.path.join('DGS', '__init__.py')) as f:
    line = f.readline()
    while not line.startswith('__version__'):
        line = f.readline()
exec(line, globals())

ext_modules = [ ]
cmdclass = { }

if USE_CYTHON:
    ext_modules += [
        Extension("DGS.cwt", [ "DGS/cwt.pyx" ],
        include_dirs=[np.get_include()]),
        Extension("DGS.sgolay", [ "DGS/sgolay.pyx" ],
        include_dirs=[np.get_include()]),
    ]
    cmdclass.update({ 'build_ext': build_ext })
else:
    ext_modules += [
        Extension("DGS.cwt", [ "DGS/cwt.c" ],
        include_dirs=[np.get_include()]),
        Extension("DGS.sgolay", [ "DGS/sgolay.c" ],
        include_dirs=[np.get_include()]),
    ]
install_requires = [
    'numpy','scipy','Pillow','matplotlib', 'cython',
]

def setupPackage():
   setup(name='pyDGS',
         version=__version__,
         description='wavelet-based digital grain size analysis',
         classifiers=[
             'Intended Audience :: Science/Research',
             'Intended Audience :: Developers',
             'License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)',
             'Programming Language :: Python',
             'Programming Language :: Python :: 2.7',
             'Programming Language :: Cython',
             'Topic :: Scientific/Engineering',
             'Topic :: Scientific/Engineering :: Physics',
         ],
         keywords='sediment',
         author='Daniel Buscombe',
         author_email='dbuscombe@usgs.gov',
         url='https://github.com/dbuscombe-usgs/pyDGS',
         download_url ='https://github.com/dbuscombe-usgs/pyDGS/archive/master.zip',
         install_requires=install_requires,
         license = "GNU GENERAL PUBLIC LICENSE v3",
         packages=['DGS'],
         cmdclass = cmdclass,
         ext_modules=ext_modules,
         platforms='OS Independent',
         include_dirs = [np.get_include()],
         package_data={'DGS': ['*.JPG','*.jpg',]}
   )

if __name__ == '__main__':
    # clean --all does not remove extensions automatically
    if 'clean' in sys.argv and '--all' in sys.argv:
        import shutil
        # delete complete build directory
        path = os.path.join(SETUP_DIRECTORY, 'build')
        try:
            shutil.rmtree(path)
        except:
            pass
        # delete all shared libs from lib directory
        path = os.path.join(SETUP_DIRECTORY, 'DGS')
        for filename in glob.glob(path + os.sep + '*.pyd'):
            try:
                os.remove(filename)
            except:
                pass
        for filename in glob.glob(path + os.sep + '*.so'):
            try:
                os.remove(filename)
            except:
                pass
    setupPackage()

