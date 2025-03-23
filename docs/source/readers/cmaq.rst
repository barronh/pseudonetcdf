.. CMAQ

CMAQ
~~~~

CMAQ uses a reader that inherits many methods from PseudoNetCDFFile, but
relies on the ioapi_base class to update meta-data.

* :class:`ioapi <PseudoNetCDF.cmaqfiles.ioapi>`: CMAQ NetCDF reader.
* :class:`cmaqomidat <PseudoNetCDF.cmaqfiles.cmaqomidat>`: OMI input file.
* :class:`griddesc <PseudoNetCDF.cmaqfiles.griddesc>`: GRIDDESC reader to make files like CMAQ that are useful as a starting place.
