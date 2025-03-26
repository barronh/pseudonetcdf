.. Readers

Formats
~~~~~~~

PseudoNetCDF helps make files act like netCDF by making readers that convert the
binary files data data and metadata, and by providing convenience functions for
the metadata. To do that, it has many `format` options. The easiest way to use
them is via :func:`pncopen`.

.. code-block:: python

    import PseudoNetCDF as pnc
    path = '/path/to/file'
    fmt = 'ioapi'  # ioapi is for CMAQ; bpch is for GEOS-Chem; CAMx has options...
    f = pnc.pncopen(path, format='?')


You can find a full list of options using `pncopen(help=True)`. The list below
gives common `format` options for specific types of files.

.. toctree::
   :maxdepth: 1
   :glob:

   readers/*


A quick summary of formats is shown below. Many formats have several aliases, and others have options.
Options might support different features (read vs read/write) or different variants of that format.
The major headers link to more information on the readers (including links to their individual
documentaton.)

+------------------------------------+-------------------------------------------------------------+
| Source                             | format text options                                         |
+====================================+=============================================================+
|            :class:`General NetCDF <PseudoNetCDF.PseudoNetCDF.core._files.netcdf>`                |
+====================================+=============================================================+
| Any NetCDF File                    | netcdf, ncf, nc, Dataset                                    |
+====================================+=============================================================+
| :doc:`CMAQ Modeling Inputs/Outputs with Coordinate (time, layer, grid) Support <readers/cmaq>`   |
+====================================+=============================================================+
| IOAPI NetCDF File                  | ioapi, cmaqfiles.ioapi                                      |
+------------------------------------+-------------------------------------------------------------+
| GRIDDESC Text File                 | griddesc, cmaqfiles.griddesc                                |
+------------------------------------+-------------------------------------------------------------+
| OMI Input                          | cmaqomidat, cmaqfiles.cmaqomidat                            |
+------------------------------------+-------------------------------------------------------------+
| BCON Input Profile                 | bcon_profile, cmaqfiles.profile.bcon_profile                |
+------------------------------------+-------------------------------------------------------------+
| ICON Input Profile                 | icon_profile, cmaqfiles.profile.icon_profile                |
+------------------------------------+-------------------------------------------------------------+
| Photolysis Table                   | jtable, cmaqfiles.jtable                                    |
+====================================+=============================================================+
| :doc:`CAMx Modeling Inputs/Outputs with Coordinate (time, layer, grid) Support <readers/camx>`   |
+====================================+=============================================================+
| Gridded Fortran File               | uamiv, camxfiles.uamiv.Read.uamiv,                          |
|                                    | camxfiles.uamiv.Memmap.uamiv                                |
+------------------------------------+-------------------------------------------------------------+
| Fine Instantaneous                 | finst, camxfiles.finst.Memmap.finst                         |
+------------------------------------+-------------------------------------------------------------+
| Integrated Reaction Rate           | irr, camxfiles.irr.Read.irr,                                |
|                                    | camxfiles.irr.Memmap.irr                                    |
+------------------------------------+-------------------------------------------------------------+
| Integrated Process Rate            | ipr, camxfiles.ipr.Read.ipr,                                |
|                                    | camxfiles.ipr.Memmap.ipr                                    |
+------------------------------------+-------------------------------------------------------------+
| Lateral Boundary Fortran           | lateral_boundary,                                           |
|                                    | camxfiles.lateral_boundary.Memmap.lateral_boundary          |
+------------------------------------+-------------------------------------------------------------+
| Wind Fortran                       | wind, camxfiles.wind.Memmap.wind,                           |
|                                    | camxfiles.wind.Read.wind                                    |
+------------------------------------+-------------------------------------------------------------+
| Vertical Diffusivity (kv)          | vertical_diffusivity,                                       |
|                                    | camxfiles.vertical_diffusivity.Memmap.vertical_diffusivity, |
|                                    | camxfiles.vertical_diffusivity.Read.vertical_diffusivity    |
+------------------------------------+-------------------------------------------------------------+
| Temperature                        |  temperature, camxfiles.temperature.Memmap.temperature,     |
|                                    | camxfiles.temperature.Read.temperature                      |
+------------------------------------+-------------------------------------------------------------+
| Point Source                       | point_source, camxfiles.point_source.Memmap.point_source,   |
|                                    | camxfiles.point_source.Read.point_source                    |
+------------------------------------+-------------------------------------------------------------+
| Landuse                            | landuse, camxfiles.landuse.Memmap.landuse                   |
+------------------------------------+-------------------------------------------------------------+
| Humidity                           | humidity, camxfiles.humidity.Memmap.humidity,               |
|                                    | camxfiles.humidity.Read.humidity                            |
+------------------------------------+-------------------------------------------------------------+
| Height/Pressure                    | height_pressure,                                            |
|                                    | camxfiles.height_pressure.Memmap.height_pressure,           |
|                                    | camxfiles.height_pressure.Read.height_pressure              |
+------------------------------------+-------------------------------------------------------------+
| Cloud/Rain                         | cloud_rain, camxfiles.cloud_rain.Memmap.cloud_rain          |
+====================================+=============================================================+
|         :doc:`GEOS-Chem Modeling Inputs and Outputs with Coordinate Support <readers/geoschem>`  |
+====================================+=============================================================+
| GEOS-Chem NetCDF                   | gcnc, geoschemfiles.gcnc                                    |
+------------------------------------+-------------------------------------------------------------+
| GEOS-Chem Binary Punch File        | bpch, bpch1, bpch2, geoschemfiles.bpch                      |
+------------------------------------+-------------------------------------------------------------+
| GEOS-Chem Flight Logs              | flightlogs, geoschemfiles.flightlogs                        |
+------------------------------------+-------------------------------------------------------------+
| GEOS Binary Output                 | geos, geoschemfiles.geos                                    |
+====================================+=============================================================+
|         Weather Research and Forecsting (WRF) Input/Ouptut With Time and Projection Support      |
+====================================+=============================================================+
| Weather Research Forecasting (WRF) | wrf, wrffiles.wrf                                           |
+====================================+=============================================================+
|         :doc:`NOAA Air Resources Laboratory HYSPLIT Model <readers/hysplit>`                     |
+====================================+=============================================================+
| Trajectory Dump Text               | arltrajdump, noaafiles.arltrajdump                          |
+------------------------------------+-------------------------------------------------------------+
| Concentration Dump Packed Bit      | arlconcdump, noaafiles.arlconcdump                          |
+------------------------------------+-------------------------------------------------------------+
| Particle Dump Dump Packed Bit      | arlpardump, noaafiles.arlpardump                            |
+------------------------------------+-------------------------------------------------------------+
| Generic Packed Bit                 | arlpackedbit, noaafiles.arlpackedbit                        |
+====================================+=============================================================+
|         :doc:`Miscellaneous Observation File Formats <readers/textfiles>`                        |
+====================================+=============================================================+
| World O3 and UV Data Centre Sonde  | woudcsonde, woudcfiles.woudcsonde                           |
+------------------------------------+-------------------------------------------------------------+
| NOAA ESRL Sonde Text               | l100, noaafiles.l100                                        |
+------------------------------------+-------------------------------------------------------------+
| AMES File Format 1001 (ICARTT)     | ffi1001, icarttfiles.ffi1001.ffi1001                        |
+------------------------------------+-------------------------------------------------------------+
| AERMOD Output                      | reader, aermodfiles.reader                                  |
+------------------------------------+-------------------------------------------------------------+
| Total Ozone (TOMS) text file       | tomsl3, toms.level3.tomsl3                                  |
+------------------------------------+-------------------------------------------------------------+
| "Raw" Format Text                  | aqsraw, epafiles.aqsraw                                     |
+------------------------------------+-------------------------------------------------------------+
| CSV text file                      | csv, textfiles.csv                                          |
+------------------------------------+-------------------------------------------------------------+
| hdf file                           | ceilometerl2, ceilometerfiles.ceilometerl2                  |
+====================================+=============================================================+
