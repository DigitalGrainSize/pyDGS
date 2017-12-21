# encoding: utf-8
"""
pyDGS - a Python framework for wavelet-based digital grain size analysis

pyDGS is an open-source project dedicated to provide a Python framework to compute estimates of grain size distribution  using the continuous wavelet transform method of Buscombe (2013) from an image of sediment where grains are clearly resolved. DOES NOT REQUIRE CALIBRATION

This program implements the algorithm of:

Buscombe, D. (2013) Transferable Wavelet Method for Grain-Size Distribution from Images of Sediment Surfaces and Thin Sections, and Other Natural Granular Patterns. Sedimentology 60, 1709-1732

https://www.danielbuscombe.com/s/Buscombe_2013_sedimentology_101111-sed12049.pdf

 Author:  Daniel Buscombe
           Northern Arizona University
           Flagstaff, AZ 86001
           daniel.buscombe@nau.edu
 Revision Dec 21, 2017
 First Revision January 18 2013

For more information visit https://github.com/dbuscombe-usgs/pyDGS

https://www.danielbuscombe.com/s/Buscombe_2013_sedimentology_101111-sed12049.pdf

Please contact:
daniel.buscombe@nau.edu

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

__version__ = '3.0.9'

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

from DGS._dgs_class_web import dgs
from DGS.test import *


