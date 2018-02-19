from unittest import TestSuite, findTestCases, TextTestRunner
test_suite = TestSuite()
def addTestCasesFromModule(module):
    test_suite.addTests(findTestCases(module))

from ..core import _files
addTestCasesFromModule(_files)

from .. import sci_var
addTestCasesFromModule(sci_var)

from .. import ArrayTransforms
addTestCasesFromModule(ArrayTransforms)

from .. import camxfiles
addTestCasesFromModule(camxfiles.wind.Memmap)
addTestCasesFromModule(camxfiles.humidity.Memmap)
addTestCasesFromModule(camxfiles.temperature.Memmap)
addTestCasesFromModule(camxfiles.vertical_diffusivity.Memmap)
addTestCasesFromModule(camxfiles.one3d.Memmap)
addTestCasesFromModule(camxfiles.height_pressure.Memmap)
addTestCasesFromModule(camxfiles.landuse.Memmap)
addTestCasesFromModule(camxfiles.cloud_rain.Memmap)
addTestCasesFromModule(camxfiles.uamiv.Memmap)
addTestCasesFromModule(camxfiles.uamiv.Write)
addTestCasesFromModule(camxfiles.point_source.Memmap)
addTestCasesFromModule(camxfiles.lateral_boundary.Memmap)
addTestCasesFromModule(camxfiles.ipr.Memmap)
addTestCasesFromModule(camxfiles.irr.Memmap)
addTestCasesFromModule(camxfiles.FortranFileUtil)
addTestCasesFromModule(camxfiles.timetuple)

from .. import net_balance
addTestCasesFromModule(net_balance)

from .. import geoschemfiles
addTestCasesFromModule(geoschemfiles._bpch)
addTestCasesFromModule(geoschemfiles._newbpch)
addTestCasesFromModule(geoschemfiles._geos)

from .. import textfiles
addTestCasesFromModule(textfiles._delimited)

from .. import icarttfiles
addTestCasesFromModule(icarttfiles.ffi1001)

from .. import _getreader
addTestCasesFromModule(_getreader)

def test_all():
	TextTestRunner(verbosity=2).run(test_suite)

if __name__=='__main__':
	test_all()

