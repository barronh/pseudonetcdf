HeadURL="$HeadURL$"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy$"
__version__ = RevisionNum

from numpy import array,where,logical_or,repeat,mean,sum,zeros
#This Package modules
from pyPA.utils.CAMxFiles import wind as reg_wind, \
					   height_pressure as reg_height_pressure, \
					   temperature as reg_temperature, \
					   cloud_rain as reg_cloud_rain, \
					   vertical_diffusivity as reg_vertical_diffusivity, \
					   humidity as reg_humidity
from pyPA.utils.sci_var import PseudoNetCDFFile, \
					PseudoNetCDFVariable, \
					PseudoNetCDFVariables, \
					PseudoNetCDFVariableConvertUnit, \
					Pseudo2NetCDF
from pyPA.utils.ArrayTransforms import CenterCAMxWind, \
							CenterTime, \
							CAMxHeightToDepth

class add_derived(PseudoNetCDFFile):
	"""add_derived provides a simple interface to add derived variables
	to a PseudoNetCDFFile interface
	
	create a new class with the following modifications:
	overwrite __childclass__ with the base class
	overwrite __addvars__ with a list of keys for variable names you 
						  intend to derive
	for each key name, create a key interface funciton (e.g. key=DEPTH, interface=__DEPTH__)
	"""
	__childclass__=None
	__addvars__={}
	def __init__(*args,**kwds):
		self=args[0]
		self.__child=self.__childclass__(*args[1:],**kwds)
		self.dimensions=self.__child.dimensions
		self.variables=PseudoNetCDFVariables(self.__variables__,self.__child.variables.keys()+self.__addvars__)
	def __variables__(self,key):
		if key in self.__addvars__:
		   return getattr(self,'__'+key+'__')()
		else:
			return self.__child.variables[key]

class time_avg_new_unit(PseudoNetCDFFile):
	__reader__=None
	def __init__(self,rffile,rows,cols,outunit={},endhour=True):
		self.__file=self.__reader__(rffile,rows,cols)
		self.dimensions={}
		self.createDimension('TSTEP',self.__file.dimensions['TSTEP']-1)
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
		return PseudoNetCDFVariableConvertUnit(self.__decorator(var,PseudoNetCDFVariable(self,k,var.typecode(),var.dimensions,values=CenterTime(var))),outunit)
	
	def __decorator(self,ovar,nvar):
		for a,v in ovar.__dict__.iteritems():
			setattr(nvar,a,v)
		return nvar

class window(PseudoNetCDFFile):
	"""
		window can convert take a rectangular prizm subset of
		dimensions and variables. metadata is unaffected.
	"""
	def __init__(self,ncffile,tslice=slice(None),kslice=slice(None),jslice=slice(None),islice=slice(None)):
		self.__file=ncffile
		self.__idx=[tslice,kslice,jslice,islice]
		self.dimensions=self.__file.dimensions.copy()
		any_non_time_key=[k for k in self.__file.variables.keys() if 'TFLAG' not in k][0]
		self.dimensions['TSTEP'],self.dimensions['LAY'],self.dimensions['ROW'],self.dimensions['COL'] \
		                    = self.__file.variables[any_non_time_key][self.__idx].shape
		self.variables=PseudoNetCDFVariables(self.__variables,self.__file.variables.keys())
		
	def __variables(self,k):
		if 'TFLAG' in k:
			return self.__file.variables[k]
		
		ov=self.__file.variables[k]
		nv=ov[self.__idx]
		Pseudo2NetCDF().addVariableProperties(nv,ov)
		return nv

class newresolution(PseudoNetCDFFile):
	"""
		newresolution can convert dimensions and variables 
		to a new resolution. metadata is unaffected.
	"""
	def __init__(self,ncffile,axis,oldres,newres,repeat_method=repeat,condense_method=sum):
		self.__axis=array(axis,ndmin=1)
		self.__oldres=array(oldres,ndmin=1)
		self.__newres=array(newres,ndmin=1)
		self.__condense=condense_method
		self.__repeat=repeat_method
		self.__file=ncffile
			
		if not logical_or(self.__oldres/self.__newres % 1 == 0,self.__newres/self.__oldres % 1 ==0).any():
			raise ValueError, "One resolution must be a factor of the other."

		self.dimensions=self.__file.dimensions.copy()
		any_non_time_key=[k for k in self.__file.variables.keys() if 'TFLAG' not in k][0]
		v=self.__file.variables[any_non_time_key]
		v=self.__method(v)
		self.dimensions['TSTEP'],self.dimensions['LAY'],self.dimensions['ROW'],self.dimensions['COL'] \
		                    =  v.shape
		self.variables=PseudoNetCDFVariables(self.__variables,self.__file.variables.keys())
	def __method(self,a):
		axis=self.__axis
		oldres=self.__oldres
		newres=self.__newres
		for i in range(axis.size):
			if oldres[i]>newres[i]:
				method=lambda a: self.__repeat(a,oldres[i]/newres[i],axis[i])
			elif newres[i]>oldres[i]:
				newshape=list(a.shape)
				newshape[axis[i]:axis[i]+1]=newshape[axis[i]]/newres[i]*oldres[i],newres[i]/oldres[i]
				method=lambda a: self.__condense(a.reshape(newshape),axis[i])
			a=method(a)
		return a

	def __variables(self,k):
		if 'TFLAG' in k and (self.__axis!=0).any():
			raise KeyError, "Tflag is off limits"
		else:
			ov=self.__file.variables[k]
			v=self.__method(ov)
			Pseudo2NetCDF().addVariableProperties(ov,v)
			return v
				
class MetaNetCDF(PseudoNetCDFFile):
	__metavars__={}
	def addMetaVariable(self,key,func):
		self.variables.addkey(key)
		self.__metavars__[key]=func

	def __init__(self,files):
		self.__files=files
		self.dimensions={}
		keys=[]
		for f in self.__files:
			for k,d in f.dimensions.iteritems():
				if d==1 and k=='LAY':
					k='SURFLAY'
				if k not in self.dimensions.keys():
					self.createDimension(k,d)
			keys.extend(f.variables.keys())
			for k in f.__dict__.keys():
				if k not in self.__dict__.keys():
					setattr(self,k,getattr(f,k))
		keys=list(set(keys))
		self.variables=PseudoNetCDFVariables(self.__variables,keys)
		self.dimensions['VAR']=len(keys)-1
	
	def __getattribute__(self,k):
		try:
			return PseudoNetCDFFile.__getattribute__(self,k)
		except AttributeError:
			for f in self.__files:
				try:
					return getattr(f,k)
				except:
					pass
			raise AttributeError, "%s not found" % k

	def __variables(self,k):
		if k in self.__metavars__.keys():
			return self.__metavars__[k](self)
		for f in self.__files:
			if k in f.variables.keys():
				v=f.variables[k]
				if k=='TFLAG':
					v=PseudoNetCDFVariable(self,'TFLAG','i',v.dimensions,values=array(v)[:,[0],:].repeat(self.dimensions['VAR'],1))
					v.long_name='TFLAG'.ljust(16)
					v.var_desc='TFLAG'.ljust(16)
					v.units='DATE-TIME'
					
				if v.shape[1]==1 and k=='LAY' and k in self.dimensions.keys():
					dims=list(v.dimensions)
					dims[1]='SURFLAY'
					v.dimensions=tuple(dims)
				return v
		raise KeyError,'%s not in any files' % k
file_master=MetaNetCDF