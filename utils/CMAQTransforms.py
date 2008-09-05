from numpy import array,ones
from pyPA.utils.MetaNetCDF import add_derived, \
                                  file_master
from pyPA.utils.ArrayTransforms import CenterTime
from pynetcdf import NetCDFFile
from pyPA.utils.ArrayTransforms import CenterCMAQWind, \
                                       CenterTime
from pyPA.utils.units import F2K
from pyPA.utils.sci_var import PseudoNetCDFFile, \
                    PseudoNetCDFVariable, \
                    PseudoIOAPIVariable, \
                    PseudoNetCDFVariables, \
                    PseudoNetCDFVariableConvertUnit
ncf=NetCDFFile
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
        return PseudoNetCDFVariableConvertUnit(self.__decorator(var,PseudoNetCDFVariable(self,k,var.typecode(),var.dimensions,values=CenterTime(var))),outunit)
    
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
            v=PseudoNetCDFVariable(self,k,'f',('TSTEP','LAY','ROW','COL'),values=preproc(var))
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
        airmols=PseudoNetCDFVariable(self,'AIRMOLS','f',('TSTEP','LAY','ROW','COL'),values=P*V/R/T)
        airmols.units='moles'
        airmols.long_name='AIRMOLS'.ljust(16)
        airmols.var_desc='AIRMOLS'.ljust(16)
        return airmols

def cmaq_pa_master(paths_and_readers,tslice=slice(None),kslice=slice(None),jslice=slice(None),islice=slice(None)):
	"""
	CMAQ PA Master presents a single interface for CMAQ PA, IRR, and 
	Instantaneous concentration files.
	
	paths_and_readers - iterable of iterables (n x 2) where each element of the
						primary iterable is an iterable containing a file path 
						and a reader for that path.  The reader is expected to
						present the Scientific.IO.NetCDF.NetCDFFile file inter-
						face.
	
	optional:
	   tslice - slice object used to window the time period from the 
				instantaneous concentration file to match the PA domain
	   kslice - same as tslice, but for layers
	   jslice - same as tslice, but for rows
	   islice -  - same as tslice, but for columns
			  
	"""
	files=[]
	for p,r in paths_and_readers:
		files.append(eval(r)(p))
	master_file=file_master(files)
	def InitLambda(x,tslice,kslice,jslice,islice):
		return lambda self: PseudoIOAPIVariable(self,x,'f',('TSTEP','LAY','ROW','COL'),values=self.variables[x][:-1,:,:,:][tslice,kslice,jslice,islice],units=self.variables[x].units)
	def FinalLambda(x,tslice,kslice,jslice,islice):
		return lambda self: PseudoIOAPIVariable(self,x,'f',('TSTEP','LAY','ROW','COL'),values=self.variables[x][1:,:,:,:][tslice,kslice,jslice,islice],units=self.variables[x].units)
	def MetLambda(x,tslice,kslice,jslice,islice):
		return lambda self: PseudoIOAPIVariable(self,x,'f',('TSTEP','LAY','ROW','COL'),values=CenterTime(self.variables[x])[tslice,kslice,jslice,islice],units=self.variables[x].units)
	for k in master_file.variables.keys():
		if '_' not in k and k!='TFLAG':
			master_file.addMetaVariable('INIT_'+k,InitLambda(k,tslice,kslice,jslice,islice))
			master_file.addMetaVariable('FCONC_'+k,FinalLambda(k,tslice,kslice,jslice,islice))
			master_file.addMetaVariable('INITIAL_'+k,InitLambda(k,tslice,kslice,jslice,islice))
			master_file.addMetaVariable('FINAL_'+k,FinalLambda(k,tslice,kslice,jslice,islice))
	master_file.addMetaVariable('CONC_AIRMOLS',lambda self: PseudoIOAPIVariable(self,'CONC_AIRMOLS','f',('TSTEP','LAY','ROW','COL'),values=CenterTime(self.variables['PRES'][:,:,:,:]/8.314472/self.variables['TA'])[tslice,kslice,jslice,islice],units='moles/m**3'))
	master_file.addMetaVariable('AIRMOLS',lambda self: PseudoIOAPIVariable(self,'AIRMOLS','f',('TSTEP','LAY','ROW','COL'),values=self.variables['CONC_AIRMOLS']*self.XCELL*self.YCELL*2*CenterTime(self.variables['ZF'][:,:,:,:]-self.variables['ZH'][:,:,:,:])[tslice,kslice,jslice,islice],units='moles'))
	master_file.addMetaVariable('INVAIRMOLS',lambda self: PseudoIOAPIVariable(self,'INVAIRMOLS','f',('TSTEP','LAY','ROW','COL'),values=1/self.variables['AIRMOLS'][:,:,:,:],units='moles'))
	master_file.addMetaVariable('DEFAULT_SHAPE',lambda self: PseudoIOAPIVariable(self,'DEFAULT_SHAPE','f',('TSTEP','LAY','ROW','COL'),values=ones(self.variables['PRES'][tslice,kslice,jslice,islice].shape,'bool'),units='on/off'))
	
	return master_file
