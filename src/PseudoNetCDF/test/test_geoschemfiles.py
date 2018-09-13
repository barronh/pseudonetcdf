__all__ = ['test_bpch', 'test_bpch1', 'test_bpch2', 'test_geos']


from PseudoNetCDF.geoschemfiles._bpch import TestMemmaps as test_bpch1
from PseudoNetCDF.geoschemfiles._newbpch import TestMemmaps as test_bpch2
from PseudoNetCDF.geoschemfiles._bpchmaster import TestMemmaps as test_bpch
from PseudoNetCDF.geoschemfiles._geos import TestMemmaps as test_geos
