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
__all__ = ['boxmodel', 'pa', 'box_model_mrg', 'box_model_conc', 'jtable',
           'bcon_profile', 'icon_profile', 'ioapi', 'ioapi_base', 'griddesc',
           'cmaqomidat']

from . import boxmodel
from . import pa
from ._jtable import jtable
from .boxmodel import box_model_conc, box_model_mrg
from .profile import bcon_profile, icon_profile
from ._ioapi import ioapi, ioapi_base
from ._griddesc import griddesc
from ._cmaqomidat import cmaqomidat
