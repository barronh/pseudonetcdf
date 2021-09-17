# PseudoNetCDF like NetCDF except for many scientific format backends

[![Build Status](https://app.travis-ci.com/barronh/pseudonetcdf.svg?branch=master)](https://app.travis-ci.com/github/barronh/pseudonetcdf)

# Overview

PseudoNetCDF provides read, plot, and sometimes write capabilities for atmospheric science data formats including:

* CAMx (www.camx.org)
* RACM2 box-model outputs
* Kinetic Pre-Processor outputs
* ICARTT Data files (ffi1001)
* CMAQ Files
* GEOS-Chem Binary Punch/NetCDF files
* and many more, for a full list of formats see `pncdump --list-formats`

Documentation on PseudoNetCDF is available at https://pseudonetcdf.readthedocs.io/

# Example code

Code below needs paths and formats filled in. The format should be chosen from "Current Formats" below.
Note that there are many methods for modifying the files. Use `help` on infile to learn more options.

```
import PseudoNetCDF as pnc
inpath = '<path-to-input>'
outpath = '<path-for-output>'
infile = pnc.pncopen(inpath, format = '<choose-format-below>')
# Print CDL representation - good for learning dimensions, variables, and properties
print(infile)
# Optionally, add dimension slicing
# infile = infile.sliceDimensions(<layer-dim-name> = 0)
# infile = infile.applyAlongDimenions(<time-dim-name> = 'mean')
# patches = infile.plot('<varkey>', plottype = 'longitude-latitude')
# patches.axes.figure.savefig('<path-for-figure>')
infile.save(outpath)
```

# Example Command Line Interfaces

pncdump - like ncdump
pncgen - reads from a file and generates a netcdf file
pncload - read in a file and open a python environment

```
$ pncdump -f <choose-format-below> <path-to-input>
$ pncgen -f <choose-format-below> <path-to-input> <path-to-output>
$ pncload -f <choose-format-below> <path-to-input>
```

# Current Formats

| Long Name | Short Name |
| ----------- | --------- |
| netcdf | netcdf |
| aermodfiles.reader | |
| cmaqfiles.ioapi | ioapi |
| cmaqfiles.jtable | jtable |
| epafiles.aqsraw | aqsraw |
| geoschemfiles.bpch | bpch |
| geoschemfiles.bpch2 | bpch2 |
| geoschemfiles.flightlogs | flightlogs |
| geoschemfiles.geos | geos |
| noaafiles.arlpackedbit | arlpackedbit |
| noaafiles.arlconcdump | arlconcdump |
| noaafiles.arlpackedbit | arlpackedbit |
| noaafiles.arlpardump | arlpardump |
| noaafiles.arltrajdump | arltrajdump |
| textfiles.csv | csv |
| cmaqfiles.profile.bcon_profile | bcon_profile |
| cmaqfiles.profile.icon_profile | icon_profile |
| icarttfiles.ffi1001.ffi1001 | ffi1001 |
| camxfiles.cloud_rain.Memmap.cloud_rain | cloud_rain |
| camxfiles.finst.Memmap.finst | finst |
| camxfiles.height_pressure.Memmap.height_pressure | height_pressure |
| camxfiles.height_pressure.Read.height_pressure | |
| camxfiles.humidity.Memmap.humidity | humidity |
| camxfiles.humidity.Read.humidity | |
| camxfiles.ipr.Memmap.ipr | ipr |
| camxfiles.ipr.Read.ipr | |
| camxfiles.irr.Memmap.irr | irr |
| camxfiles.irr.Read.irr | |
| camxfiles.landuse.Memmap.landuse | landuse |
| camxfiles.lateral_boundary.Memmap.lateral_boundary | lateral_boundary |
| camxfiles.point_source.Memmap.point_source | point_source |
| camxfiles.point_source.Read.point_source | |
| camxfiles.temperature.Memmap.temperature | temperature |
| camxfiles.temperature.Read.temperature | |
| camxfiles.uamiv.Memmap.uamiv | uamiv |
| camxfiles.uamiv.Read.uamiv | |
| camxfiles.uamiv.Transforms.osat | osat |
| camxfiles.vertical_diffusivity.Memmap.vertical_diffusivity | vertical_diffusivity |
| camxfiles.vertical_diffusivity.Read.vertical_diffusivity | |
| camxfiles.wind.Memmap.wind | wind |
| camxfiles.wind.Read.wind | |

# More information

Lots more available at our [wiki ](http://github.com/barronh/pseudonetcdf/wiki)

Try our:
  * [Install Instructions](http://github.com/barronh/pseudonetcdf/wiki/Install-Instructions)
  * [CAMx Tutorials](http://github.com/barronh/pseudonetcdf/wiki/CAMx-Tutorials)
  * [GEOS-Chem Tutorials](http://github.com/barronh/pseudonetcdf/wiki/GC-Tutorials)
  * [Recipes](Recipes)


Quick tour:
 * [Install](http://github.com/barronh/pseudonetcdf/wiki/Install-Instructions.md):
  * `pip install https://github.com/barronh/pseudonetcdf/archive/v3.1.0.zip` for the most stable version or 
  * `pip install http://github.com/barronh/pseudonetcdf/archive/master.zip` for the latest.
 * Download example icartt file: e.g., [HOx from INTEX-NA](http://www-air.larc.nasa.gov/cgi-bin/enzFile?c16141B08DF7F1ACFBAD5C83F9313E20C792f7075622d6169722f4152435441532f4443385f41495243524146542f4252554e452e57494c4c49414d2f484f785f4443385f32303038303632365f52312e696374)
  * `curl -L ftp://ftp-air.larc.nasa.gov/pub/INTEXA/DC8_AIRCRAFT/BRUNE.WILLIAM/HOX_DC8_20040626_R0.ict`
 * Dump an icartt file in CDL: `pncdump -f ffi1001 HOX_DC8_20040626_R0.ict`
 * Create a netcdf from an icartt file: `pncgen -f ffi1001 HOX_DC8_20040626_R0.ict HOX_DC8_20040626_R0.nc`
