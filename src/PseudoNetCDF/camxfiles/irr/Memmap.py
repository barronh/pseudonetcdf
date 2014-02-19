__all__=['irr']
__doc__ = """
.. _Memmap
:mod:`Memmap` -- irr Memmap interface
============================================

.. module:: Memmap
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` memory map for CAMx
              irr files.  See PseudoNetCDF.sci_var.PseudoNetCDFFile 
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
from warnings import warn
from collections import defaultdict
#Site-Packages
from numpy import zeros,array,where,memmap,newaxis,dtype,nan

#This Package modules
from PseudoNetCDF.camxfiles.timetuple import timediff,timeadd,timerange
from PseudoNetCDF.camxfiles.FortranFileUtil import OpenRecordFile,Int2Asc
from PseudoNetCDF.sci_var import PseudoNetCDFFile, PseudoNetCDFVariable, PseudoNetCDFVariables
from PseudoNetCDF.ArrayTransforms import ConvertCAMxTime

#for use in identifying uncaught nan
listnan=struct.unpack('>f','\xff\xc0\x00\x00')[0]
checkarray=zeros((1,),'f')
checkarray[0]=listnan
array_nan=checkarray[0]

def _isMine(path):
    testf = OpenRecordFile(path)
    if testf.record_size == 80:
        testf.next()
        if testf.record_size == 16:
            testf.next()
            if testf.record_size == 28:
                return True
    return False

class irr(PseudoNetCDFFile):
    """
    irr provides a PseudoNetCDF interface for CAMx
    irr files.  Where possible, the inteface follows
    IOAPI conventions (see www.baronams.com).
    
    ex:
        >>> irr_path = 'camx_irr.bin'
        >>> rows,cols = 65,83
        >>> irrfile = irr(irr_path,rows,cols)
        >>> irrfile.variables.keys()
        ['TFLAG', 'RXN_01', 'RXN_02', 'RXN_03', ...]
        >>> v = irrfile.variables['RXN_01']
        >>> tflag = irrfile.variables['TFLAG']
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
        >>> irrfile.dimensions
        {'TSTEP': 25, 'LAY': 28, 'ROW': 65, 'COL': 83}
    """
    
    id_fmt="if5i"
    data_fmt="f"
    @classmethod
    def isMine(self, path):
        return _isMine(path)
    def __init__(self,rf,multi=False):
        """
        Initialization included reading the header and learning
        about the format.
        
        see __readheader and __gettimestep() for more info
        """
        self.__rffile=OpenRecordFile(rf)
        self.__readheader()
        irr_record_type=dtype(
                    dict(names=(['SPAD','DATE', 'TIME', 'PAGRID', 'NEST', 'I', 'J', 'K']+
                                [ 'RXN_%02d' % i for i in range(1,self.nrxns+1)]+
                                ['EPAD']),
                         formats=(['>i', '>i', '>f', '>i', '>i', '>i', '>i', '>i']+ 
                                 [ '>f' for i in range(1,self.nrxns+1)]+
                                 ['>i'])))
    
        varkeys=[i for i in irr_record_type.names[8:-1]]+['TFLAG']
        self.groups = defaultdict(PseudoNetCDFFile)
        padatatype = []
        pavarkeys = []
        self.createDimension('TSTEP', self.time_step_count)
        self.createDimension('VAR',  len(varkeys)-1)
        self.createDimension('DATE-TIME', 2)
        for di, domain in enumerate(self._padomains):
            dk ='PA%02d' % (di + 1)
            prefix = dk + '_'
            pavarkeys.extend([prefix + k for k in varkeys])
            grp = self.groups[dk]
            for propk, propv in domain.iteritems():
                setattr(grp, propk, propv)
            grp.createDimension('TSTEP', self.time_step_count)
            grp.createDimension('VAR',  len(varkeys)-1)
            grp.createDimension('DATE-TIME', 2)
            grp.createDimension('COL', domain['iend']-domain['istart']+1)
            grp.createDimension('ROW', domain['jend']-domain['jstart']+1)
            grp.createDimension('LAY', domain['tlay']-domain['blay']+1)
            padatatype.append((dk, irr_record_type, (len(grp.dimensions['ROW']), len(grp.dimensions['COL']), len(grp.dimensions['LAY']))))
            if len(self._padomains) == 1:
                self.createDimension('COL', domain['iend']-domain['istart']+1)
                self.createDimension('ROW', domain['jend']-domain['jstart']+1)
                self.createDimension('LAY', domain['tlay']-domain['blay']+1)
                for propk, propv in domain.iteritems():
                    setattr(grp, propk, propv)
            exec("""def varget(k):
                return self._irr__variables('%s', k)""" % dk, dict(self = self), locals())
            if len(self._padomains) == 1:
                self.variables = PseudoNetCDFVariables(varget,varkeys)
            else:
                grp.variables = PseudoNetCDFVariables(varget,varkeys)
        self.__memmaps = memmap(self.__rffile.infile.name, dtype(padatatype), 'r', self.data_start_byte)

    def __del__(self):
        try:
            self.__memmaps.close()
            del self.__memmaps
        except:
            pass
                
    
    def __decorator(self,name,pncfv):
        decor=lambda k: dict(units='ppm/hr', var_desc=k.ljust(16), long_name=k.ljust(16))
        for k,v in decor(name).iteritems():
            setattr(pncfv,k,v)        
        return pncfv
        
    def __variables(self,pk, rxn):
        if rxn=='TFLAG':
            return ConvertCAMxTime(self.__memmaps[pk][:, 0,0,0]['DATE'],self.__memmaps[pk][:,0,0,0]['TIME'],len(self.groups[pk].dimensions['VAR']))
        return self.__decorator(rxn,PseudoNetCDFVariable(self,rxn,'f',('TSTEP','LAY','ROW','COL'),values=self.__memmaps[pk][rxn].swapaxes(1, 3).swapaxes(2, 3)))
        
    def __readheader(self):
        """
        __readheader reads the header section of the ipr file
        it initializes each header field (see CAMx Users Manual for a list)
        as properties of the ipr class
        """
        assert(self.__rffile.record_size == 80)
        self.runmessage=self.__rffile.read("80s")
        self.start_date,self.start_time,self.end_date,self.end_time=self.__rffile.read("ifif")
        self.time_step=100.
        self.time_step_count=len([i for i in self.timerange()])
        self._grids=[]
        for grid in range(self.__rffile.read("i")[-1]):
            self._grids.append(
                            dict(
                                zip(
                                    ['orgx','orgy','ncol','nrow','xsize','ysize','iutm'], 
                                    self.__rffile.read("iiiiiii")
                                    )
                                )
                            )
        
        self._padomains=[]
        for padomain in range(self.__rffile.read("i")[-1]):
            self._padomains.append(
                                dict(
                                    zip(
                                        ['grid','istart','iend','jstart','jend','blay','tlay'],
                                        self.__rffile.read("iiiiiii")
                                        )
                                    )
                                )
        self.nrxns=self.__rffile.read('i')[-1]
        
        self.data_start_byte=self.__rffile.record_start
        self.record_fmt=self.id_fmt + str(self.nrxns) + self.data_fmt
        self.record_size=self.__rffile.record_size
        self.padded_size=self.record_size+8
        domain=self._padomains[0]
        self.records_per_time=(domain['iend']-domain['istart']+1)*(domain['jend']-domain['jstart']+1)*(domain['tlay']-domain['blay']+1)
        self.time_data_block=self.padded_size*self.records_per_time
        self.time_step=100.

    def timerange(self):
        return timerange((self.start_date,self.start_time+self.time_step),timeadd((self.end_date,self.end_time),(0,self.time_step)),self.time_step)

class TestMemmap(unittest.TestCase):
    def runTest(self):
        pass
    def setUp(self):
        pass
        
    def testIRR(self):
        warn('Test not implemented')
       
if __name__ == '__main__':
    unittest.main()
