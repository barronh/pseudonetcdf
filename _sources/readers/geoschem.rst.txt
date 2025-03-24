.. GEOS-Chem

GEOS-Chem
~~~~~~~~~

GEOS-Chem used the Binary Punch (bpch) file for many years, but now relies on NetCDF.
PseudoNetCDF provides readers and convenience functions for both.

* :class:`gcnc <PseudoNetCDF.geoschemfiles.gcnc>`: Newer GEOS-Chem NetCDF files
* :class:`bpch <PseudoNetCDF.geoschemfiles.bpch>`: Best Binary Punch reader (chooses bpch1 or bpch2)
* :class:`bpch1 <PseudoNetCDF.geoschemfiles.bpch1>`: Binary Punch reader (load to mem)
* :class:`bpch2 <PseudoNetCDF.geoschemfiles.bpch2>`: Binary Punch reader (memory mapped)
* :class:`geos <PseudoNetCDF.geoschemfiles.geos>`: Raw GEOS met files.
* :class:`flightlogs <PseudoNetCDF.PseudoNetCDF.geoschemfiles._planelog.flightlogs>`: Legacy flight log output.
