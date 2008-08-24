from numpy import array
from pyPA.utils.MetaNetCDF import add_derived, \
                                  file_master
from pynetcdf import NetCDFFile
from pyPA.utils.ArrayTransforms import CenterCMAQWind, \
                                       CenterTime
from pyPA.utils.units import F2K
from pyPA.utils.sci_var import PseudoNetCDFFile, \
                    PseudoNetCDFVariable, \
                    PseudoNetCDFVariables, \
                    PseudoNetCDFVariableConvertUnit

#==================================================================
#                                                             time_avg_new_unit 
class time_avg_new_unit(PseudoNetCDFFile):
#    __reader__=None
    def __init__(self,rffile,outunit={},endhour=False):
        self.__file=NetCDFFile(rffile)
        self.dimensions={}
        self.createDimension('TSTEP',self.__file.variables[self.__file.variables.keys()[0]].shape[0]-1)
        self.createDimension('LAY',self.__file.dimensions['LAY'])
        self.createDimension('ROW',self.__file.dimensions['ROW'])
        self.createDimension('COL',self.__file.dimensions['COL'])
        self.createDimension('VAR',self.__file.dimensions['VAR'])
        self.createDimension('DATE-TIME',self.__file.dimensions.get('DATE-TIME',2))
        self.__outunit=outunit
        self.variables=PseudoNetCDFVariables(self.__variables,self.__file.variables.keys())
        self.__timeslice={True:slice(1,None),False:slice(None,-1)}[endhour]
        v=self.createVariable('TFLAG','i',('TSTEP','VAR','DATE-TIME'),keep=True)
        v.assignValue(self.__file.variables['TFLAG'][self.__timeslice])
        v.long_name='Time flag'
        v.units='DATE-TIME'

    def __variables(self,k):
        outunit=self.__outunit.get(k,None)
        var=self.__file.variables[k]
        if outunit==None:
            outunit=var.units
        return PseudoNetCDFVariableConvertUnit(self.__decorator(var,PseudoNetCDFVariable(self,k,var.typecode(),var.dimensions,CenterTime(var))),outunit)
    
    def __decorator(self,ovar,nvar):
        for a,v in ovar.__dict__.iteritems():
            setattr(nvar,a,v)
        return nvar

#==================================================================
#                                                             NetCDFFile_center_time 
# class NetCDFFile_center_time(time_avg_new_unit):
#     __reader__=NetCDFFile

#==================================================================
#                                                             NetCDFFile_center_time 
class wind_center_time_cell(PseudoNetCDFFile):
    """
    CMAQ Files
    """
    def __init__(self,rffile,outunit='m/s',endhour=False):
        self.__windfile=NetCDFFile(rffile)
        self.dimensions={}
        self.createDimension('TSTEP',self.__windfile.variables['UWIND'].shape[0]-1)
        self.createDimension('LAY',self.__windfile.dimensions['LAY'])
        self.createDimension('ROW',self.__windfile.dimensions['ROW']-1)
        self.createDimension('COL',self.__windfile.dimensions['COL']-1)
        self.createDimension('VAR',self.__windfile.dimensions['VAR'])
        self.createDimension('DATE-TIME',self.__windfile.dimensions.get('DATE-TIME',2))
        self.__outunit=outunit
        self.variables=PseudoNetCDFVariables(self.__variables,['VWIND','UWIND','TFLAG'])
        self.__timeslice={True:slice(1,None),False:slice(None,-1)}[endhour]
    def __variables(self,k):
        self.__add_variables()
        return self.variables[k]
        
    def __add_variables(self):
        v=self.createVariable('TFLAG','i',('TSTEP','VAR','DATE-TIME'),keep=True)
        v.assignValue(self.__windfile.variables['TFLAG'][self.__timeslice])
        v.long_name='Time flag'
        v.units='DATE-TIME'
        for k in ['UWIND','VWIND']:
            preproc=CenterCMAQWind   
            var=self.__windfile.variables[k]
            v=PseudoNetCDFVariable(self,k,'f',('TSTEP','LAY','ROW','COL'),preproc(var))
            v.units=var.units
            v.long_name=k.ljust(16)
            v.var_desc=(k+' at center').ljust(16)
            self.variables[k]=PseudoNetCDFVariableConvertUnit(v,self.__outunit)
    
    

#==================================================================

def pypass_cmaq_met_master(metcro2d_path,metcro3d_path,metdot3d_path,rows,cols,endhour=True):
     windf=wind_center_time_cell(metdot3d_path,outunit='km/h',endhour=endhour)
     dim2f=time_avg_new_unit(metcro2d_path,outunit={'PBL':'km', 'TEMPG':'deg_F'},endhour=endhour)
     dim3f=time_avg_new_unit(metcro3d_path,outunit={'TA':'deg_F','ZF':'km'},endhour=endhour)
     return file_master([windf,dim2f,dim3f])
     
def pypass_cmaq_emiss_master(emiss_path,rows,cols,endhour=True):
     emiss = time_avg_new_unit(emiss_path,outunit={'CO':'moles/h', 'NO':'moles/h', 'NO2':'moles/h', 'ALD2':'moles/h', 'ETH':'moles/h','ETOH':'moles/h', 'MEOH':'moles/h','FORM':'moles/h', 'ISOP':'moles/h', 'OLE':'moles/h', 'PAR':'moles/h', 'TOL':'moles/h', 'XYL':'moles/h'},endhour=endhour)
     return file_master([emiss])

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