PseudoNetCDF provides read, plot, and sometimes write capabilities for atmospheric science data formats including:

* CAMx (www.camx.org)
* RACM2 box-model outputs
* Kinetic Pre-Processor outputs
* ICARTT Data files (ffi1001)
* CMAQ Files
* GEOS-Chem Binary Punch/NetCDF files

for a full list of formats see `pncdump --help` and look in the --formats help section.

Try our:
  * [Install Instructions](http://github.com/barronh/pseudonetcdf/wiki/Install-Instructions)
  * [CAMx Tutorials](http://github.com/barronh/pseudonetcdf/wiki/CAMx-Tutorials)
  * [GEOS-Chem Tutorials](http://github.com/barronh/pseudonetcdf/wiki/GC-Tutorials)
  * [Recipes](Recipes)


Quick tour:
 * [Install](http://github.com/barronh/pseudonetcdf/wiki/Install-Instructions.md): `pip install git+git://github.com/barronh/pseudonetcdf.git`
 * Download example icartt file: e.g., [HOx from INTEX-NA](ftp://ftp-air.larc.nasa.gov/pub/INTEXA/DC8_AIRCRAFT/BRUNE.WILLIAM/HOX_DC8_20040626_R0.ict)
 * Dump an icartt file in CDL: `pncdump -f ffi1001 HOX_DC8_20040626_R0.ict`
 * Create a netcdf from an icartt file: `pncgen -f ffi1001 HOX_DC8_20040626_R0.ict HOX_DC8_20040626_R0.nc`
