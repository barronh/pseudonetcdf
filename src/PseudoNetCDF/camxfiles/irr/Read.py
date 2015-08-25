__all__=['irr']
__doc__ = """
.. _Read
:mod:`Read` -- irr Read interface
============================================

.. module:: Read
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` random access read for CAMx
              irr files.  See PseudoNetCDF.sci_var.PseudoNetCDFFile 
              for interface details
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""
HeadURL="$HeadURL: http://dawes.sph.unc.edu:8080/uncaqmlsvn/pyPA/utils/trunk/CAMxRead.py $"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy: svnbarronh $"
__version__ = RevisionNum

#Distribution packages
from datetime import datetime, timedelta
from types import GeneratorType
import unittest
import struct,sys,os,operator
from warnings import warn
from tempfile import TemporaryFile as tempfile
from math import ceil
import os,sys

#Site-Packages
from numpy import zeros,array,where,memmap,newaxis,dtype, fromfile, arange

#This Package modules
from PseudoNetCDF.camxfiles.util import sliceit
from PseudoNetCDF.sci_var import PseudoNetCDFFile, PseudoNetCDFVariable, PseudoNetCDFVariables


#for use in identifying uncaught nan
listnan=struct.unpack('>f','\xff\xc0\x00\x00')[0]
checkarray=zeros((1,),'f')
checkarray[0]=listnan
array_nan=checkarray[0]

class irr(PseudoNetCDFFile):
    """
    irr provides a PseudoNetCDF interface for CAMx
    irr files.  Where possible, the inteface follows
    IOAPI conventions (see www.baronams.com).
    
    ex:
        >>> irr_path = 'camx_irr.bin'
        >>> irrfile = irr(irr_path)
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
    
    id_fmt="ifiiiii"
    data_fmt="f"
    def __init__(self,rf,units='ppm/hr',conva=None, nvarcache = None):
        """
        Initialization included reading the header and learning
        about the format.
        
        see __readheader and __gettimestep() for more info
        """
        rffile = self._rffile = open(rf)
        rffile.seek(0,2)
        if rffile.tell()<2147483648L:
            warn("For greater speed on files <2GB use ipr_memmap")
        rffile.seek(0,0)
        self.__readheader()
        self.__gettimestep()
        self.units=units
        #__conv is a conversion array that comes from ipr
        #if units is not ppm/hr, conv must be provided
        self.__conv=conva
        if self.units!='ppm/hr' and self.__conv==None:
            raise ValueError, "When units are provided, a conversion array dim(t,z,x,y) must also be provided"
        varkeys=['RXN_%02d' % i for i in range(1,self.NRXNS+1)]

        domain=self.padomains[0]
        tdim = self.createDimension('TSTEP', int(self.NSTEPS))
        tdim.setunlimited(True)
        self.createDimension('COL', int(domain['iend']-domain['istart']+1))
        self.createDimension('ROW', int(domain['jend']-domain['jstart']+1))
        self.createDimension('LAY', int(domain['tlay']-domain['blay']+1))
        self.createDimension('DATE-TIME', 2)
        self.createDimension('VAR', int(self.NRXNS))
        self.variables=PseudoNetCDFVariables(self.__var_get,varkeys)
        self._nvarcache = nvarcache or int(ceil(self.NRXNS / 4.))
        tflag = self.createVariable('TFLAG', 'i', ('TSTEP', 'VAR', 'DATE-TIME'))
        tflag.units = '<YYYYJJJ, HHMMSS>'
        tflag.var_desc = tflag.long_name = 'TFLAG'.ljust(16)
        tflag[:] = array(self.timerange(), dtype = 'i').reshape(self.NSTEPS, 1, 2)

    def __var_get(self,key):
        rxni = int(key.split('_')[1])
        self.loadVars(rxni, self._nvarcache)
        tflag = self.createVariable('TFLAG', 'i', ('TSTEP', 'VAR', 'DATE-TIME'))
        tflag.units = '<YYYYJJJ, HHMMSS>'
        tflag.var_desc = tflag.long_name = 'TFLAG'.ljust(16)
        tflag[:] = array(self.timerange(), dtype = 'i').reshape(self.NSTEPS, 1, 2)
        return self.variables[key]

    def __readheader(self):
        """
        __readheader reads the header section of the ipr file
        it initializes each header field (see CAMx Users Manual for a list)
        as properties of the ipr class
        """
        rffile = self._rffile
        self.runmessage=fromfile(rffile, dtype = dtype([('SPAD', '>i'), ('NOTE', 'S80'), ('EPAD', '>i')]), count = 1)[0]['NOTE']
        self.SDATE, self.STIME, self.EDATE, self.ETIME = fromfile(rffile, dtype = dtype([('SPAD', '>i'), ('SDATE', '>i'), ('STIME', '>f'), ('EDATE', '>i'), ('ETIME', '>f'), ('EPAD', '>i')]), count = 1)[0].tolist()[1:-1]
        
        self.grids=[]
        for grid in range(fromfile(rffile, dtype = dtype([('SPAD', '>i'), ('NGRIDS', '>i'), ('EPAD', '>i')]), count = 1)[0]['NGRIDS']):
            self.grids.append(
                            dict(
                                zip(
                                    ['orgx','orgy','ncol','nrow','xsize','ysize','iutm'], 
                                    fromfile(rffile, dtype = dtype([('SPAD', '>i'), ('orgx', '>f'), ('orgy', '>f'), ('ncol', '>i'), ('nrow', '>i'), ('xsize', '>f'), ('ysize', '>f'), ('iutm', '>i'), ('EPAD', '>i')]), count = 1)[0].tolist()[1:-1]
                                    )
                                )
                            )
        
        self.padomains=[]
        for padomain in range(fromfile(rffile, dtype = dtype([('SPAD', '>i'), ('NPADOMAINS', '>i'), ('EPAD', '>i')]), count = 1)[0]['NPADOMAINS']):
            self.padomains.append(
                                dict(
                                    zip(
                                        ['grid','istart','iend','jstart','jend','blay','tlay'],
                                        fromfile(rffile, dtype = dtype([('SPAD', '>i'), ('grid', '>i'), ('istart', '>i'), ('iend', '>i'), ('jstart', '>i'), ('jend', '>i'), ('blay', '>i'), ('tlay', '>i'), ('EPAD', '>i')]), count = 1)[0].tolist()[1:-1]
                                        )
                                    )
                                )
        self.NRXNS = fromfile(rffile, dtype = dtype([('SPAD', '>i'), ('NRXNS', '>i'), ('EPAD', '>i')]), count = 1)[0]['NRXNS']
        
        self._data_start_byte=rffile.tell()
        
        self._record_dtype = dtype([('SPAD', '>i'), ('DATE', '>i'), ('TIME', '>f'), ('PAGRID', '>i'), ('NEST', '>i'), ('I', '>i'), ('J', '>i'), ('K', '>i'), ('IRRS', '>f', self.NRXNS), ('EPAD', '>i')])
        self._padded_size=self._record_dtype.itemsize
        
    def __gettimestep(self):
        """
        Header information provides start and end date, but does not
        indicate the increment between.  This routine reads the first
        and second date/time and initializes variables indicating the
        timestep length and the anticipated number.
        """
        self._activedomain=self.padomains[0]
        self._rffile.seek(
                        self._data_start_byte+(
                                    self.__jrecords(0,self.padomains[0]['jend'])*
                                    self._padded_size
                                    ), 0
                        )
        date,time= fromfile(self._rffile, dtype = self._record_dtype, count = 1)[0].tolist()[1:3]
        tstart = datetime.strptime('%05dT%04d' % (self.SDATE, self.STIME), '%y%jT%H%M')
        tone = datetime.strptime('%05dT%04d' % (date, time), '%y%jT%H%M')
        tstep = tone - tstart
        self.time_step = self.TSTEP = int((datetime.strptime('0000', '%H%M') + tstep).strftime('%H%M'))
        tend = datetime.strptime('%05dT%04d' % (self.EDATE, self.ETIME), '%y%jT%H%M')
        tdiff = tend - tstart
        multiple = (tdiff.days * 24 * 3600. + tdiff.seconds) / (tstep.days * 24 * 3600. + tstep.seconds)
        self.NSTEPS = self.time_step_count = int(multiple)
        assert(multiple == int(multiple))

    def __gridrecords(self,pagrid):
        """
        routine returns the number of records to increment from the
        data start byte to find the pagrid
        """
        ntime=self.__timerecords(pagrid,(self.EDATE,int(self.ETIME+self.TSTEP)))
        return ntime
        
    def __timerecords(self,pagrid,(d,t)):
        """
        routine returns the number of records to increment from the
        data start byte to find the first time
        """
        nsteps=self.timerange().index((d,t))
        nj=self.__jrecords(pagrid,self.padomains[pagrid]['jend']+1)
        return nsteps*nj
        
    def __irecords(self,pagrid,i):
        """
        routine returns the number of records to increment from the
        data start byte to find the first icell
        """
        ni=i-self._activedomain['istart']
        nk=self.__krecords(pagrid,self._activedomain['tlay']+1)
        return ni*nk
        
    def __jrecords(self,pagrid,j):
        """
        routine returns the number of records to increment from the
        data start byte to find the first jcell
        """
        nj=j-self._activedomain['jstart']
        ni=self.__irecords(pagrid,self._activedomain['iend']+1)
        return nj*ni
        
    def __krecords(self,pagrid,k):
        """
        routine returns the number of records to increment from the
        data start byte to find the first kcell
        """
        return k-self._activedomain['blay']

    def __recordposition(self,pagrid,date,time,i,j,k):
        """ 
        routine uses pagridrecords, timerecords,irecords,
        jrecords, and krecords multiplied by the fortran padded size
        to return the byte position of the specified record
        
        pagrid - integer
        date - integer
        time - float
        i - integer
        j - integer
        k - integer
        """
        records=0
        for pag in range(pagrid):
            records+=__gridrecords(pag)
        records+=self.__timerecords(pagrid,(date,time))
        records+=self.__jrecords(pagrid,j)
        records+=self.__irecords(pagrid,i)
        records+=self.__krecords(pagrid,k)
        return records*self._padded_size+self._data_start_byte
        
    def _seek(self,pagrid=1,date=None,time=None,i=1,j=1,k=1):
        """
        Move file cursor to beginning of specified record
        see __recordposition for a definition of variables
        """
        if date==None:
            date=self.SDATE
        if time==None:
            time=self.STIME+self.TSTEP
        self._activedomain=self.padomains[pagrid]
        self._rffile.seek(self.__recordposition(pagrid,date,time,i,j,k))
        
    def loadVars(self,start, n, pagrid=0):
        domain=self.padomains[pagrid]
        istart=domain['istart']
        iend=domain['iend']
        jstart=domain['jstart']
        jend=domain['jend']
        kstart=domain['blay']
        kend=domain['tlay']
        variables = self.variables
        nk = kend + 1 - kstart
        nj = jend + 1 - jstart
        ni = iend + 1 - istart
        nrec = nk * ni * nj
        temp = zeros((nrec, self.NRXNS), 'f')
        shape = (self.NSTEPS,) + tuple(eval('map(len, (LAY, ROW, COL))', None, self.dimensions))
        variables.clear()
        end = min(start + n, self.NRXNS + 1)
        start = max(1, end - n)
        for rxn in range(start, end):
            key = 'RXN_%02d' % rxn
            variables[key] = PseudoNetCDFVariable(self, key, 'f', ('TSTEP', 'LAY', 'ROW', 'COL'), values = zeros(shape, 'f'), units = 'ppm/hr', var_desk = key.ljust(16), long_name = key.ljust(16))

        self._seek(pagrid = 0, i = istart, j = jstart, k = kstart)
        for ti, (d,t) in enumerate(self.timerange()):
            record = fromfile(self._rffile, dtype = self._record_dtype, count = nrec)
            temp[:] = record['IRRS']
            date = record['DATE']
            time = record['TIME']
            id = record['I']
            jd = record['J']
            kd = record['K']
            assert((id == arange(istart, iend + 1)[None, :, None].repeat(nk, 2).repeat(nj, 0).ravel()).all())
            assert((kd == arange(kstart, kend + 1)[None, None, :].repeat(ni, 1).repeat(nj, 0).ravel()).all())
            assert(((jd == arange(jstart, jend + 1).repeat(ni * nk, 0))).all())
            assert((date == d).all() and (time == t).all())
            for rxn in range(start, end):
                variables['RXN_%02d' % rxn][ti, :, :, :] = temp[:, rxn-1].reshape(nj, ni, nk).swapaxes(1,2).swapaxes(0, 1)
            
    def timerange(self):
        tstart = datetime.strptime('%05dT%04d' % (self.SDATE, self.STIME), '%y%jT%H%M')
        tdiff = datetime.strptime('%04d' % self.TSTEP, '%H%M') - datetime.strptime('0000', '%H%M')
        dates = [tstart + (tdiff * i) for i in range(1, self.NSTEPS+1)]
        return [(int(d.strftime('%y%j')), float(d.strftime('%H%M'))) for d in dates]

class TestRead(unittest.TestCase):
    def runTest(self):
        pass
    def setUp(self):
        pass
        
    def testIRR(self):
        emissfile=irr('../../../../testdata/ei/camx_cb4_ei_lo.20000825.hgb8h.base1b.psito2n2.hgbpa_04km')
        self.assert_(1==2)
       
if __name__ == '__main__':
    unittest.main()
