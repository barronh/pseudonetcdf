__doc__ = """
.. _one3d
:mod:`one3d` -- one3d File Interfaces
=====================================

.. module:: one3d
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` memory map and random access read 
   based file interfaces for CAMx one variable 3d files.
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""
__all__=['Memmap','Read','Write']

import Memmap
import Read
import Write
