__all__=['vertical_diffusivity']
__doc__ = """
.. _Read
:mod:`Read` -- vertical/diffusivity Read interface
============================================

.. module:: Read
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` random access read for CAMx
              vertical/diffusivity files.  See PseudoNetCDF.sci_var.PseudoNetCDFFile 
              for interface details
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""
HeadURL="$HeadURL: http://dawes.sph.unc.edu:8080/uncaqmlsvn/pyPA/utils/trunk/CAMxRead.py $"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy: svnbarronh $"
__version__ = RevisionNum

#Distribution packages
from types import GeneratorType
import unittest
import struct,sys,os,operator
from warnings import warn
from tempfile import TemporaryFile as tempfile
import os,sys

#Site-Packages
from numpy import zeros,array,where,memmap,newaxis,dtype

#This Package modules
from PseudoNetCDF.camxfiles.timetuple import timediff,timeadd,timerange
from PseudoNetCDF.camxfiles.FortranFileUtil import OpenRecordFile,read_into,Int2Asc,Asc2Int
from PseudoNetCDF.sci_var import PseudoNetCDFFile, PseudoNetCDFVariable, PseudoNetCDFVariables
from PseudoNetCDF.camxfiles.one3d.Read import one3d

#for use in identifying uncaught nan
listnan=struct.unpack('>f','\xff\xc0\x00\x00')[0]
checkarray=zeros((1,),'f')
checkarray[0]=listnan
array_nan=checkarray[0]

class vertical_diffusivity(one3d):
    """
    vertical_diffusivity provides a PseudoNetCDF interface for CAMx
    vertical_diffusivity files.  Where possible, the inteface follows
    IOAPI conventions (see www.baronams.com).
    
    ex:
        >>> vd_path = 'camx_vd.bin'
        >>> rows,cols = 65,83
        >>> vdfile = vertical_diffusivity(vd_path,rows,cols)
        >>> vdfile.variables.keys()
        ['TFLAG', 'KV']
        >>> tflag = vdfile.variables['TFLAG']
        >>> tflag.dimensions
        ('TSTEP', 'VAR', 'DATE-TIME')
        >>> tflag[0,0,:]
        array([2005185,       0])
        >>> tflag[-1,0,:]
        array([2005185,  240000])
        >>> v = vdfile.variables['KV']
        >>> v.dimensions
        ('TSTEP', 'LAY', 'ROW', 'COL')
        >>> v.shape
        (25, 28, 65, 83)
        >>> vdfile.dimensions
        {'TSTEP': 25, 'LAY': 28, 'ROW': 65, 'COL': 83}
    """
    
    var_name='KV'
    units='m**2/s'


class TestRead(unittest.TestCase):
    def runTest(self):
        pass
    def setUp(self):
        pass
        
    def testKV(self):
        vdfile=vertical_diffusivity('../../../../testdata/met/camx_kv.20000825.hgbpa_04km.TCEQuh1_eta.v43.tke',65,83)
        self.assert_((vdfile.variables['KV'].mean(0).mean(1).mean(1)==array([  13.65080357,   34.39198303,   68.02783966,   95.5898819 ,
          109.25765991,  112.92014313,  108.32209778,   97.25794983,
          84.1328125 ,   65.92033386,   46.97774506,   25.8343792 ,
          9.80327034,    2.89653206,    1.26993668,    1.12098336,
          1.13557184,    1.13372564,    1.19559622,    1.1675849 ,
          1.18877947,    1.18713808,    1.02371764,    1.02544105,
          1.21638143,    1.34624374,    1.03213251,    1.        ],dtype='f')).all())

if __name__ == '__main__':
    unittest.main()
