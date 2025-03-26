.. CAMx

CAMx
~~~~

CAMx has many readers for inputs, but the uamiv reader covers all the outputs.
If you're using CAMx with NetCDF outputs, you can use the CMAQ reader.


* New NetCDF Inpput/Output Files:

    * :class:`ioapi <PseudoNetCDF.cmaqfiles.ioapi>`: NetCDF IOAPI metadata reader, or
    * :class:`netcdf <PseudoNetCDF.PseudoNetCDF.core._files.netcdf>`: Standard NetCDF.


* Fortran Formatted Outputs:

    * :class:`uamiv <PseudoNetCDF.camxfiles.Memmaps.uamiv>`: CAMx (and UAMIV) format reads
       * average and instantaneous output
       * gridded emission and initial condition input
    * :class:`finst <PseudoNetCDF.camxfiles.Memmaps.finst>`: Finescale instantaneous outputs
    * :class:`ipr <PseudoNetCDF.camxfiles.Memmaps.ipr>`: Integrated Process Rate (ipr) output file
    * :class:`irr <PseudoNetCDF.camxfiles.Memmaps.irr>`: Integrated Reaction Rate (ipr) output file

* Fortran Binary Inputs:
    * :class:`cloud_rain <PseudoNetCDF.camxfiles.Memmaps.cloud_rain>`: Cloud/rain input file
    * :class:`height_pressure <PseudoNetCDF.camxfiles.Memmaps.height_pressure>`: Height/Pressure (zp) file
    * :class:`humidity <PseudoNetCDF.camxfiles.Memmaps.humidity>`: Humidity input file
    * :class:`landuse <PseudoNetCDF.camxfiles.Memmaps.landuse>`: Landuse input file
    * :class:`lateral_boundary <PseudoNetCDF.camxfiles.Memmaps.lateral_boundary>`: Lateral boundary condition file
    * :class:`point_source <PseudoNetCDF.camxfiles.Memmaps.point_source>`: Point source input file
    * :class:`temperature <PseudoNetCDF.camxfiles.Memmaps.temperature>`: Temperature file
    * :class:`vertical_diffusivity <PseudoNetCDF.camxfiles.Memmaps.vertical_diffusivity>`: Vertical diffusivity (KV) file.
    * :class:`wind <PseudoNetCDF.PseudoNetCDF.camxfiles.wind.Memmap.wind>`: Wind input file.
