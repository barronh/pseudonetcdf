In the sciences, file formats abound and it is easy to spend more effort reading/writing data than analyzing it. The goal of PseudoNetCDF is to provide a single abstract interface for many data formats. The project is inspired by [NetCDF](http://www.unidata.ucar.edu/software/netcdf/) PseudoNetCDF abstracts many file formats using the data/object model from NetCDF.  The NetCDF library has a strong object model, a long history, and a long list of meta-data conventions.  Rather than designing our own interface, we rely on the power and simplicity of the Scientific.IO.NetCDF python interface for the NetCDF library.  PseudoNetCDF provides this same data interface for the following data formats from CAMx (www.camx.org), RACM2 box-model outputs, Kinetic Pre-Processor outputs, CMAQ box-model outputs, ICARTT formatted data, and more.

Try our:
  * [Install Instructions](InstallInstructions.md)
  * [Basic Tutorials](BasicTutorials.md)
  * [CAMx Tutorials](CAMxTutorials.md)
  * [Recipes](Recipes.md)


Quick tour:
  * Install: pip install git+https://code.google.com/p/pseudonetcdf/
  * Download example icartt file: e.g., [HOx fromm INTEX-NA](ftp://ftp-air.larc.nasa.gov/pub/INTEXA/DC8_AIRCRAFT/BRUNE.WILLIAM/HOX_DC8_20040626_R0.ict)
  * Dump an icartt file in CDL: python -m PseudoNetCDF.pncdump -f ffi1001 HOX\_DC8\_20040626\_R0.ict
  * Create a netcdf from an icartt file: python -m PseudoNetCDF.pncgen -f ffi1001 HOX\_DC8\_20040626\_R0.ict HOX\_DC8\_20040626\_R0.nc