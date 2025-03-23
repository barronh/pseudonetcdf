.. PseudoNetCDF documentation master file, created by
   sphinx-quickstart on Tue Sep  4 21:46:43 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PseudoNetCDF User's Guide
=========================

About
-----

The key value of PseudoNetCDF is netcdf-like access to complex data formats
in the air quality field.

PseudoNetCDF was created to provide netcdf-like (see netCDF4) access air
quality data written in other formats. The interface was designed based on
Scientific.io.netcdf and netcdf4-python. Then it grew to provide meta-data
aware spatial and temporal processing, some inspired by xarray and pandas.

* :doc:`installing`
* :doc:`quick`: Short examples.
* :doc:`examples`: Short real-world applications.
* :doc:`readers`: Details about available formats.
* :doc:`core`: Some high-level description of what makes PseuodNetCDF work.

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Table of Contents

   self
   quick
   examples
   readers
   installing
   core
   api/PseudoNetCDF
   issues
