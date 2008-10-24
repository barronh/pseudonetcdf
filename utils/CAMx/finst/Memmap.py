HeadURL="$HeadURL: http://dawes.sph.unc.edu:8080/uncaqmlsvn/pyPA/utils/trunk/CAMxMemmap.py $"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy: svnbarronh $"
__version__ = RevisionNum

__all__=['uamiv']
#Distribution packages
import unittest
import struct

#Site-Packages
from numpy import zeros,array,where,memmap,newaxis,dtype,nan,fromfile

#This Package modules
from pyPA.utils.timetuple import timediff,timeadd
from pyPA.utils.FortranFileUtil import OpenRecordFile,Int2Asc
from pyPA.utils.sci_var import PseudoNetCDFFile, PseudoNetCDFVariable, PseudoNetCDFVariables
from pyPA.utils.ArrayTransforms import ConvertCAMxTime

#for use in identifying uncaught nan
listnan=struct.unpack('>f','\xff\xc0\x00\x00')[0]
checkarray=zeros((1,),'f')
checkarray[0]=listnan
array_nan=checkarray[0]

class finst(PseudoNetCDFFile):
	"""
	This class is intended to provide an interface to the
	nested initial conditions outputs/innputs of the 
	CAMx model
	"""
	
	__messagefmt='>240S'
	__nestspcfmt=dtype(dict(names=['nnest','nspec'],formats=['>i','>i']))
	__nestparmsfmt=dtype(dict(names=['SPAD','ibeg','jbeg','iend','jend','mesh','ione','nx','ny','nz','iparent','ilevel','EPAD'],formats=['>i']*13))
	__time_hdr_fmt=dtype(dict(names=['time','date'],formats=['>f','>i']))
	__spc_fmt=dtype(">10S")
	__ione=1
	__idum=0
	__rdum=0.
	__buffersize=4
	def __init__(self,rf):
		"""
		Initialization included reading the header and learning
		about the format.
		
		see __readheader and __gettimestep() for more info
		"""
		self.__rffile=rf

		# Establish dimensions
		self.dimensions={'DATE-TIME': 2 }

		self.__readheader()
		
		# Add IOAPI metavariables
		nlays=self.NLAYS=self.dimensions['LAY']
		nrows=self.NROWS=self.dimensions['ROW']
		ncols=self.NCOLS=self.dimensions['COL']
		nvars=self.NVARS=self.dimensions['VAR']
		nsteps=self.NSTEPS=self.dimensions['TSTEP']
		setattr(self,'VAR-LIST',"".join([i.ljust(16) for i in self.__var_names__] + ['TFLAG'.ljust(16)]))
		self.GDTYP=2

		# Create variables
		self.variables=PseudoNetCDFVariables(self.__variables,self.__var_names__+['TFLAG'])

		# Initialize time maps
		date=self.__memmap__['DATE']
		time=self.__memmap__['TIME']
		self.variables['TFLAG']=ConvertCAMxTime(date,time,self.NVARS)
		self.SDATE,self.STIME=self.variables['TFLAG'][0,0,:]
		

	def __checkfilelen(self):
		f=file(self.__rffile,'rb')
		f.seek(0,2)
		flen=f.tell()
		f.close()
		return flen

	def __readheader(self):
		start=0
		end=0
		#start+=self.__buffersize
		#end=self.__messagefmt.itemsize/4
		#self.MESSAGE=self.__memmap__[start:end].view(self.__messagefmt)
		f=file(self.__rffile)
		f.seek(92)
		self.NNEST,self.NSPEC=fromfile(f,'>i',2)
		f.seek(8,1)

		self.SPECIES=fromfile(f,self.__spc_fmt,self.NSPEC)
		
		offset=f.tell()+4
		f.close()
		del f
		self.__var_names__=[i.strip() for i in self.SPECIES.tolist()]

		self.__memmap__=memmap(self.__rffile,'>f','r',offset=offset)

		start=0 #skip end and start buffer
		end=start+self.__nestparmsfmt.itemsize/4*self.NNEST		
		self.NEST_PARMS=self.__memmap__[start:end].view(self.__nestparmsfmt)

		date_time_block_size=4
		spc_1_lay_block_size=self.NEST_PARMS['nx']*self.NEST_PARMS['ny']+2
		grid_block_sizes=(self.NSPEC*self.NEST_PARMS['nz']*spc_1_lay_block_size)
		grid_fmt=[]
		for i in range(self.NNEST):
			grid_fmt.append('>%df' % int(grid_block_sizes[i]))
		data_block_size=date_time_block_size+grid_block_sizes.sum()

		ntimes=float(self.__memmap__.size-end)/data_block_size
		if int(ntimes)!=ntimes:
			raise ValueError, "Not an even number of times"
		ntimes=int(ntimes)

		self.createDimension('TSTEP',ntimes)
		self.createDimension('VAR',self.NSPEC)
		self.chooseGrid(0)
		
		start=end
		end=start+ntimes*data_block_size
		if end!=self.__memmap__.size:
			raise ValueError, "Dimensions do not match file"
		self.__memmap__=self.__memmap__[start:].view(dtype(dict(names=['SPAD','TIME','DATE','EPAD']+['grid%d' % i for i in range(self.NNEST)],formats=['>i','>f','>i','>i']+grid_fmt)))
		
	def chooseGrid(self,ngrid):
		self.createDimension('LAY',self.NEST_PARMS[ngrid]['nz'])
		self.createDimension('COL',self.NEST_PARMS[ngrid]['nx'])
		self.createDimension('ROW',self.NEST_PARMS[ngrid]['ny'])
		self.CURRENT_GRID='grid%d' % ngrid

	def __decorator(self,name,pncfv):
		if name=='EMISSIONS ':
			decor=lambda spc: dict(units='umol/hr', var_desc=spc.ljust(16), long_name=spc.ljust(16))
		else:
			decor=lambda spc: dict(units='ppm', var_desc=spc.ljust(16), long_name=spc.ljust(16))

		for k,v in decor(name).iteritems():
			setattr(pncfv,k,v)
		return pncfv
		
	def __variables(self,k):
		spc_index=self.__var_names__.index(k)
		dimensions=('TSTEP','LAY','ROW','COL')
		ntimes=self.dimensions['TSTEP']
		nx=self.dimensions['COL']
		ny=self.dimensions['ROW']
		nz=self.dimensions['LAY']
		outvals=self.__memmap__[self.CURRENT_GRID].view('>f').reshape(ntimes,self.NSPEC,nz,ny*nx+2)[:,:,:,1:-1].reshape(ntimes,self.NSPEC,nz,ny,nx)[:,spc_index,:,:,:]
		return PseudoNetCDFVariable(self,k,'f',dimensions,values=outvals,units='ppm')

class TestMemmap(unittest.TestCase):
    def runTest(self):
        pass
    def setUp(self):
        pass
    def testGE(self):
        import pyPA.testcase
        emissfile=uamiv(pyPA.testcase.CAMxAreaEmissions)
        emissfile.variables['TFLAG']
        v=emissfile.variables['NO']
        self.assert_((emissfile.variables['NO'].mean(1).mean(1).mean(1)==array([  52.05988312,   51.58646774,   51.28796387,   55.63090134,
         63.95315933,  105.3456192 ,  158.26776123,  152.04057312,
         147.32403564,  154.80661011,  164.03274536,  171.88658142,
         174.36567688,  180.03359985,  173.81938171,  180.50257874,
         178.56637573,  161.35736084,  110.38669586,   97.90225983,
         89.08138275,   81.10474396,   73.36611938,   58.82622528],dtype='f')).all())

    def testAvg(self):
        import pyPA.testcase
        emissfile=uamiv(pyPA.testcase.CAMxAverage)
        emissfile.variables['TFLAG']
        v=emissfile.variables['NO']
        self.assert_((v.mean(1).mean(1).mean(1)==array([  9.44490694e-06,   2.17493564e-07,   6.08432686e-07,   9.48155161e-07,
         1.15099192e-05,   1.02132122e-04,   2.57815613e-04,   3.35910037e-04,
         3.17813188e-04,   2.51695659e-04,   1.85225872e-04,   1.40698961e-04,
         1.16110547e-04,   1.04519037e-04,   1.00367179e-04,   9.81789271e-05,
         8.98482831e-05,   6.31201983e-05,   2.18762198e-05,   1.78832056e-06,
         1.20556749e-07,   1.57714638e-07,   1.82648236e-07,   2.02759026e-07],dtype='f')).all())

    def testInst(self):
        from warnings import warn
        warn("Instantaneous file test not implemented")
       
if __name__ == '__main__':
    unittest.main()
