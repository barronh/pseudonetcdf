__all__=['humidity']
__doc__ = """
.. _Memmap
:mod:`Memmap` -- humidity Memmap interface
============================================

.. module:: Memmap
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` memory map for CAMx
              humidity files.  See PseudoNetCDF.sci_var.PseudoNetCDFFile 
              for interface details
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""
HeadURL="$HeadURL: http://dawes.sph.unc.edu:8080/uncaqmlsvn/pyPA/utils/trunk/CAMxMemmap.py $"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy: svnbarronh $"
__version__ = RevisionNum

#Distribution packages
import unittest
import struct

#Site-Packages
from numpy import zeros,array,where,memmap,newaxis,dtype,nan

#This Package modules
from PseudoNetCDF.camxfiles.timetuple import timediff,timeadd
from PseudoNetCDF.camxfiles.FortranFileUtil import OpenRecordFile,Int2Asc
from PseudoNetCDF.camxfiles.one3d.Memmap import one3d
from PseudoNetCDF.ArrayTransforms import ConvertCAMxTime

#for use in identifying uncaught nan
listnan=struct.unpack('>f','\xff\xc0\x00\x00')[0]
checkarray=zeros((1,),'f')
checkarray[0]=listnan
array_nan=checkarray[0]
    
class humidity(one3d):
    """
    humidity provides a PseudoNetCDF interface for CAMx
    humidity files.  Where possible, the inteface follows
    IOAPI conventions (see www.baronams.com).
    
    ex:
        >>> humidity_path = 'camx_humidity.bin'
        >>> rows,cols = 65,83
        >>> humidityfile = humidity(humidity_path,rows,cols)
        >>> humidityfile.variables.keys()
        ['TFLAG', 'HUM']
        >>> v = humidityfile.variables['HUM']
        >>> tflag = humidityfile.variables['TFLAG']
        >>> tflag.dimensions
        ('TSTEP', 'VAR', 'DATE-TIME')
        >>> tflag[0,0,:]
        array([2005185,       0])
        >>> tflag[-1,0,:]
        array([2005185,  240000])
        >>> v.dimensions
        ('TSTEP', 'LAY', 'ROW', 'COL')
        >>> v.shape
        (25, 28, 65, 83)
        >>> humidityfile.dimensions
        {'TSTEP': 25, 'LAY': 28, 'ROW': 65, 'COL': 83}
    """
    
    var_name='HUM'
    units='ppm'

class TestMemmap(unittest.TestCase):
    def runTest(self):
        pass
    def setUp(self):
        pass

    def testHUM(self):
        import PseudoNetCDF.testcase
        humfile=humidity(PseudoNetCDF.testcase.camxfiles_paths['humidity'],4,5)
        humfile.variables['TFLAG']
        self.assert_((humfile.variables['HUM']==array([3.38721924e+01, 3.40657959e+01, 3.41392822e+01, 3.42358398e+01, 3.42543945e+01, 3.38868408e+01, 3.40622559e+01, 3.42358398e+01, 3.44768066e+01, 3.46112061e+01, 3.37558594e+01, 3.39323730e+01, 3.42663574e+01, 3.46854248e+01, 3.48144531e+01, 3.39472656e+01, 3.41900635e+01, 3.46160889e+01, 3.48209229e+01, 3.47874756e+01, 8.18692688e+02, 8.31423950e+02, 8.35911926e+02, 8.41317688e+02, 8.42154602e+02, 8.19450806e+02, 8.31168518e+02, 8.42139709e+02, 8.57428589e+02, 8.65724854e+02, 8.10228210e+02, 8.21361023e+02, 8.42930908e+02, 8.71891541e+02, 8.80794861e+02, 8.22559692e+02, 8.37773499e+02, 8.67419128e+02, 8.81912903e+02, 8.79253235e+02, 6.78652344e+01, 6.82532959e+01, 6.84020996e+01, 6.85950928e+01, 6.86304932e+01, 6.78945312e+01, 6.82465820e+01, 6.85941162e+01, 6.90783691e+01, 6.93474121e+01, 6.76313477e+01, 6.79859619e+01, 6.86558838e+01, 6.94960938e+01, 6.97552490e+01, 6.80159912e+01, 6.85028076e+01, 6.93570557e+01, 6.97674561e+01, 6.97009277e+01, 3.38721924e+01, 3.40657959e+01, 3.41392822e+01, 3.42358398e+01, 3.42543945e+01, 3.38868408e+01, 3.40622559e+01, 3.42358398e+01, 3.44768066e+01, 3.46112061e+01, 3.37558594e+01, 3.39323730e+01, 3.42663574e+01, 3.46854248e+01, 3.48144531e+01, 3.39472656e+01, 3.41900635e+01, 3.46160889e+01, 3.48209229e+01, 3.47874756e+01, 8.18909973e+02, 8.31444763e+02, 8.36043457e+02, 8.42020325e+02, 8.42962769e+02, 8.19678955e+02, 8.31006470e+02, 8.42196716e+02, 8.58089050e+02, 8.66791565e+02, 8.10734680e+02, 8.21652527e+02, 8.43171753e+02, 8.72346680e+02, 8.81473267e+02, 8.23317078e+02, 8.38317505e+02, 8.67773499e+02, 8.82287720e+02, 8.79729065e+02, 6.78652344e+01, 6.82532959e+01, 6.84020996e+01, 6.85950928e+01, 6.86304932e+01, 6.78945312e+01, 6.82465820e+01, 6.85941162e+01, 6.90783691e+01, 6.93474121e+01, 6.76313477e+01, 6.79859619e+01, 6.86558838e+01, 6.94960938e+01, 6.97552490e+01, 6.80159912e+01, 6.85028076e+01, 6.93570557e+01, 6.97674561e+01, 6.97009277e+01],dtype='f').reshape(2,3,4,5)).all())

               
if __name__ == '__main__':
    unittest.main()
