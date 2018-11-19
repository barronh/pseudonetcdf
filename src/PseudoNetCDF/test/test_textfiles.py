__all__ = ['test_csv']
import os
import unittest
from PseudoNetCDF.textfiles._delimited import TestCsv as test_csv
from PseudoNetCDF.testcase import all_paths


class test_toms(unittest.TestCase):
    def setUp(self):
        self.path = all_paths['tomsl3']

    def testToms(self):
        import PseudoNetCDF as pnc
        f = pnc.pncopen(self.path, format='cdtoms')
        assert(str(f) == '')
