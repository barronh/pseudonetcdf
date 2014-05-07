from unittest import TestSuite, findTestCases, TextTestRunner
test_suite = TestSuite()
def addTestCasesFromModule(module):
    test_suite.addTests(findTestCases(module))

import sci_var
addTestCasesFromModule(sci_var)

import ArrayTransforms
addTestCasesFromModule(ArrayTransforms)

import camxfiles
addTestCasesFromModule(camxfiles.wind.Memmap)
addTestCasesFromModule(camxfiles.humidity.Memmap)
addTestCasesFromModule(camxfiles.temperature.Memmap)
addTestCasesFromModule(camxfiles.vertical_diffusivity.Memmap)
addTestCasesFromModule(camxfiles.height_pressure.Memmap)
addTestCasesFromModule(camxfiles.landuse.Memmap)
addTestCasesFromModule(camxfiles.cloud_rain.Memmap)
addTestCasesFromModule(camxfiles.uamiv.Memmap)
addTestCasesFromModule(camxfiles.point_source.Memmap)
addTestCasesFromModule(camxfiles.lateral_boundary.Memmap)
addTestCasesFromModule(camxfiles.ipr.Memmap)
addTestCasesFromModule(camxfiles.irr.Memmap)
addTestCasesFromModule(camxfiles.FortranFileUtil)
addTestCasesFromModule(camxfiles.timetuple)

import net_balance
addTestCasesFromModule(net_balance)

import geoschemfiles
addTestCasesFromModule(geoschemfiles._bpch)
addTestCasesFromModule(geoschemfiles._geos)


def test():
	TextTestRunner(verbosity=2).run(test_suite)

if __name__=='__main__':
	test()

