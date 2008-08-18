from numpy import array
from pyPA.utils.MetaNetCDF import add_derived,file_master
from pynetcdf import NetCDFFile
from pyPA.utils.units import F2K
from pyPA.utils.sci_var import PseudoNetCDFVariables,PseudoNetCDFFile,PseudoNetCDFVariable

class MetaMetPlusAirMols(add_derived):
    __childclass__=lambda *args: NetCDFFile(*args[1:])
    __addvars__=['AIRMOLS','AREA']
    def __AREA__(self):
        #self.XSIZE*self.YSIZE
        return 4000**2
        
    def __AIRMOLS__(self):
        V=array(self.variables['DEPTH'])*array(self.variables['AREA']) #m**3
        T=F2K(array(self.variables['AIRTEMP'])) #K
        P=array(self.variables['PRES'])*100. #PA
        
        R=8.314472 #m**3 x Pa x K**-1 x mol**-1
        airmols=PseudoNetCDFVariable(self,'AIRMOLS','f',('TSTEP','LAY','ROW','COL'),P*V/R/T)
        airmols.units='moles'
        airmols.long_name='AIRMOLS'.ljust(16)
        airmols.var_desc='AIRMOLS'.ljust(16)
        return airmols
        

class CMAQMetaPA(PseudoNetCDFFile):
	"""CMAQMetaPA provides a single interface to 
	multiple CMAQ PA and concentration files
	"""
	def __init__(self,papaths,concpath,kslice=slice(None),jslice=slice(None),islice=slice(None)):
		"""
		papaths - an iterable of file paths for PA or IPR files
		concpath - a single file path to an instantaneous concentration
		kslice - a vertical slice of the concentration file whose domain 
		         is a superset of the PA domain (optional)
		jslice - same as kslice, but for rows
		islice - same as kslice, but for cols
		"""
		files=[]
		PAProcess=[]
		PASpecies=[]
		self.__kslice=kslice
		self.__jslice=jslice
		self.__islice=islice
		for papath in papaths:
			files.append(NetCDFFile(papath,'r+'))
			PAProcess+=[k.split('_')[0] for k in files[-1].variables.keys() if k!='TFLAG']
			PASpecies+=[k.split('_')[1] for k in files[-1].variables.keys() if k!='TFLAG']
		PAProcess=list(set(PAProcess))
		PASpecies=list(set(PASpecies))
		files.append(NetCDFFile(concpath,'r+'))
		self.__child=file_master(files)
		self.dimensions=self.__child.dimensions
		self.variables=PseudoNetCDFVariables(self.__variables__,self.__child.variables.keys()+['INIT_'+spc for spc in PASpecies]+['FINAL_'+spc for spc in PASpecies])

	def __variables__(self,key):
		if key[:5] in ['INITI','INIT_']:
		   return self.__INITCONC__(key)
		if key[:6] in ['FCONC_','FINAL_']:
		   return self.__FINALCONC__(key)
		else:
			return self.__child.variables[key]
	def __INITCONC__(self,key):
		var=self.variables[key[5:]]
		newvar=PseudoNetCDFVariable(self,key,var.typecode(),('TSTEP','LAY','ROW','COL'),var[:-1,self.__kslice,self.__jslice,self.__islice])
		for k,v in var.__dict__.iteritems():
			setattr(newvar,k,v)
		return newvar

	def __FINALCONC__(self,key):
		var=self.variables[key[6:]]
		newvar=PseudoNetCDFVariable(self,key,var.typecode(),('TSTEP','LAY','ROW','COL'),var[1:,self.__kslice,self.__jslice,self.__islice])
		for k,v in var.__dict__.iteritems():
			setattr(newvar,k,v)
		return newvar