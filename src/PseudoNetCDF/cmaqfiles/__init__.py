__doc__ = """
.. _cmaqfiles
:mod:`cmaqfile` -- CMAQ File Interfaces
===================================

.. module:: CMAQ
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF`  random access read 
              based file interfaces for all CMAQ related files.  See
              PseudoNetCDF.sci_var.PseudoNetCDFFile for interface details.
              Each camx file level module (i.e. uamiv, cloud_rain, etc) can
              be called as a script for ncdump style output.
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""
__all__=['boxmodel', 'pa', 'box_model_mrg', 'box_model_conc', 'bcon_profile', 'icon_profile']

import boxmodel
import pa
from _jtable import jtable
from boxmodel import *
from profile import *