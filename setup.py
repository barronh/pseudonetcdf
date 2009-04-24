from distutils.core import setup
from setuptools import setup
import os
import sys
from warnings import warn
netcdf_pkgs = [('pynetcdf', 'NetCDFFile'), ('netCDF3', 'Dataset'), \
               ('netCDF4', 'Dataset'), ('Scientific.IO.NetCDF', 'NetCDFFile'), \
               ('pupynere', 'NetCDFFile')]
for pkg, reader in netcdf_pkgs:
    try:
        NetCDFFile = getattr(__import__(pkg, fromlist = [reader]),reader)
        print >> file(os.path.join('src', 'PseudoNetCDF', 'netcdf.py'),'wb'), """
__all__ = ['NetCDFFile']
__doc__ = \"\"\"
.. _netcdf
:mod:`netcdf` -- netcdf import point
====================================

.. module:: netcdf
   :platform: Unix, Windows
   :synopsis: Povides a single import point for a package.  If
              a user has one of many netcdf interfaces, this module
              selects it and provides it.
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
\"\"\"
from %s import %s as NetCDFFile
""" % (pkg,reader)
        break
    except ImportError, e:
        warn(e.message)
else:
    raise ImportError, "Did not find a NetCDFFile object"

setup(name = 'PseudoNetCDF',
      version = '1.0',
      author = 'Barron Henderson',
      author_email = 'barronh@gmail.com',
      maintainer = 'Barron Henderson',
      maintainer_email = 'barronh@gmail.com',
      packages = ['PseudoNetCDF', 'PseudoNetCDF/camxfiles', 'PseudoNetCDF/camxfiles/cloud_rain', 'PseudoNetCDF/camxfiles/finst', 'PseudoNetCDF/camxfiles/height_pressure', 'PseudoNetCDF/camxfiles/humidity', 'PseudoNetCDF/camxfiles/ipr', 'PseudoNetCDF/camxfiles/irr', 'PseudoNetCDF/camxfiles/landuse', 'PseudoNetCDF/camxfiles/one3d', 'PseudoNetCDF/camxfiles/point_source', 'PseudoNetCDF/camxfiles/temperature', 'PseudoNetCDF/camxfiles/uamiv', 'PseudoNetCDF/camxfiles/vertical_diffusivity', 'PseudoNetCDF/camxfiles/wind', 'PseudoNetCDF/cmaqfiles', 'PseudoNetCDF/racmfiles],
      package_dir = {'': 'src'},
      requires = [pkg, 'numpy (>=0.9)', 'yaml']
      )
