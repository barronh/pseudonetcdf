__all__ = ['NetCDFFile']
__doc__ = """
.. _netcdf
:mod:`netcdf` -- netcdf import point
====================================

.. module:: netcdf
   :platform: Unix, Windows
   :synopsis: Povides a single import point for a package.  If
              a user has one of many netcdf interfaces, this module
              selects it and provides it.
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""
from warnings import warn
netcdf_pkgs = [('pynetcdf', 'NetCDFFile'), ('netCDF3', 'Dataset'), \
               ('netCDF4', 'Dataset'), ('Scientific.IO.NetCDF', 'NetCDFFile'), \
               ('pupynere', 'NetCDFFile')]
for pkg, reader in netcdf_pkgs:
    try:
        NetCDFFile = getattr(__import__(pkg, fromlist = [reader]),reader)
        break
    except ImportError, e:
        #warn(e.message)
        pass
else:
    raise ImportError, "Did not find a NetCDFFile object"
