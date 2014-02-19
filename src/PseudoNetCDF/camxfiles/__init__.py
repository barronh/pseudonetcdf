__doc__ = """
.. _CAMx
:mod:`CAMx` -- CAMx File Interfaces
===================================

.. module:: CAMx
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` memory map and random access read 
              based file interfaces for all CAMx related files.  See
              PseudoNetCDF.sci_var.PseudoNetCDFFile for interface details.
              Each camx file level module (i.e. uamiv, cloud_rain, etc) can
              be called as a script for ncdump style output.  Example:
              python -m PseudoNetCDF.camxfiles.uamiv
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""
__all__=['Memmaps',
    'Readers',
    'cloud_rain',
    'height_pressure',
    'humidity',
    'ipr',
    'irr',
    'landuse',
    'lateral_boundary',
    'one3d',
    'point_source',
    'temperature',
    'uamiv',
    'vertical_diffusivity',
    'wind',
    'FortranFileUtil',
    'util',
    'timetuple']

import cloud_rain
import height_pressure
import humidity
import ipr
import irr
import landuse
import one3d
import point_source
import temperature
import uamiv
import vertical_diffusivity
import wind
import FortranFileUtil
import util
import timetuple
import lateral_boundary
import Readers
import Memmaps
