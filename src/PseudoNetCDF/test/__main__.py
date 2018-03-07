from unittest import TestSuite, findTestCases, TextTestRunner
from ..core import _files
from .. import sci_var
from .. import ArrayTransforms
from .. import camxfiles
from .. import net_balance
from .. import geoschemfiles
from .. import textfiles
from .. import icarttfiles
from .. import _getreader


test_suite = TestSuite()


def addTestCasesFromModule(module):
    test_suite.addTests(findTestCases(module))


# Core file tests
addTestCasesFromModule(_files)

# sci_var tests
addTestCasesFromModule(sci_var)

# ArrayTransforms ests
addTestCasesFromModule(ArrayTransforms)

# CAMx tests
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

# Net balance test
addTestCasesFromModule(net_balance)

# GEOS-Chem tests
addTestCasesFromModule(geoschemfiles._bpch)
addTestCasesFromModule(geoschemfiles._newbpch)
addTestCasesFromModule(geoschemfiles._geos)

# Text file tests
addTestCasesFromModule(textfiles._delimited)

# ICARTT Files
addTestCasesFromModule(icarttfiles.ffi1001)

# Reader tests
addTestCasesFromModule(_getreader)


def test_all():
    TextTestRunner(verbosity=2).run(test_suite)


if __name__ == '__main__':
    test_all()
