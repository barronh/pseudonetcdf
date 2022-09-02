import unittest
import PseudoNetCDF as pnc


class PNCVerstionTest(unittest.TestCase):
    def testVersionString(self):
        assert isinstance(pnc.__version__, str)
