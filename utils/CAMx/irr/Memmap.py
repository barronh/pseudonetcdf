HeadURL="$HeadURL: http://dawes.sph.unc.edu:8080/uncaqmlsvn/pyPA/utils/trunk/CAMxMemmap.py $"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy: svnbarronh $"
__version__ = RevisionNum

__all__=['irr']
#Distribution packages
import unittest
import struct
from warnings import warn

#Site-Packages
from numpy import zeros,array,where,memmap,newaxis,dtype,nan

#This Package modules
from pyPA.utils.timetuple import timediff,timeadd,timerange
from pyPA.utils.FortranFileUtil import OpenRecordFile,Int2Asc
from pyPA.utils.sci_var import PseudoNetCDFFile, PseudoNetCDFVariable, PseudoNetCDFVariables
from pyPA.utils.ArrayTransforms import ConvertCAMxTime

#for use in identifying uncaught nan
listnan=struct.unpack('>f','\xff\xc0\x00\x00')[0]
checkarray=zeros((1,),'f')
checkarray[0]=listnan
array_nan=checkarray[0]

class irr(PseudoNetCDFFile):
    id_fmt="if5i"
    data_fmt="f"
    def __init__(self,rf,multi=False):
        """
        Initialization included reading the header and learning
        about the format.
        
        see __readheader and __gettimestep() for more info
        """
        self.__rffile=OpenRecordFile(rf)
        self.__readheader()
        self.irr_record_type=dtype(
                        dict(names=(['SPAD','DATE', 'TIME', 'PAGRID', 'NEST', 'I', 'J', 'K']+
                                    [ 'RXN_%02d' % i for i in range(1,self.nrxns+1)]+
                                    ['EPAD']),
                                formats=(['>i', '>i', '>f', '>i', '>i', '>i', '>i', '>i']+ 
                                        [ '>f' for i in range(1,self.nrxns+1)]+
                                        ['>i'])))
    
        varkeys=[i for i in self.irr_record_type.names[8:-1]]+['TFLAG']

        domain=self.padomains[0]
        self.dimensions=dict(TSTEP=self.time_step_count,COL=domain['iend']-domain['istart']+1,ROW=domain['jend']-domain['jstart']+1,LAY=domain['tlay']-domain['blay']+1,VAR=len(varkeys)-1)
        self.variables=PseudoNetCDFVariables(self.__variables,varkeys)
        self.__memmaps=memmap(self.__rffile.infile.name,self.irr_record_type,'r',self.data_start_byte).reshape(self.dimensions['TSTEP'],self.dimensions['ROW'],self.dimensions['COL'],self.dimensions['LAY']).swapaxes(1,2).swapaxes(1,3)

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
        
    def __variables(self,rxn):
    	if rxn=='TFLAG':
    		return ConvertCAMxTime(self.__memmaps[:,0,0,0]['DATE'],self.__memmaps[:,0,0,0]['TIME'],self.dimensions['VAR'])
        return self.__decorator(rxn,PseudoNetCDFVariable(self,rxn,'f',('TSTEP','LAY','ROW','COL'),self.__memmaps[:,:,:,:][rxn]))

    def __readheader(self):
        """
        __readheader reads the header section of the ipr file
        it initializes each header field (see CAMx Users Manual for a list)
        as properties of the ipr class
        """
        self.runmessage=self.__rffile.read("80s")
        self.start_date,self.start_time,self.end_date,self.end_time=self.__rffile.read("ifif")
        self.time_step=100.
        self.time_step_count=len([i for i in self.timerange()])
        self.grids=[]
        for grid in range(self.__rffile.read("i")[-1]):
            self.grids.append(
                            dict(
                                zip(
                                    ['orgx','orgy','ncol','nrow','xsize','ysize','iutm'], 
                                    self.__rffile.read("iiiiiii")
                                    )
                                )
                            )
        
        self.padomains=[]
        for padomain in range(self.__rffile.read("i")[-1]):
            self.padomains.append(
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
        domain=self.padomains[0]
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
