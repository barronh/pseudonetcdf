__all__=['osat']
from numpy import array

from pyPA.utils.sci_var import PseudoNetCDFFile, PseudoNetCDFVariables, PseudoNetCDFVariable
from pyPA.utils.CAMxFiles import uamiv

class osat(PseudoNetCDFFile):
	__delim='_'
	def __init__(self,rffile,groups,regions):
		self.__child=uamiv(rffile)
		self.__groupsbyId=dict([('%03d' % (i+1),group) for i,group in enumerate(groups)])
		self.__groupsbyId['000']='XXX'
		self.__regionsbyId=dict([('%03d' % (i+1),region) for i,region in enumerate(regions)])
		self.__regionsbyId['IC']='IC'
		self.__regionsbyId['BC']='BC'
		self.__regionsbyNm=dict([(v,k) for k,v in self.__regionsbyId.iteritems()])
		self.__groupsbyNm=dict([(v,k) for k,v in self.__groupsbyId.iteritems()])
		self.dimensions=self.__child.dimensions.copy()
		spc_keys=[k for k in self.__child.variables.keys() if k[:3] in ['NOX','VOC','O3N','O3V']]
		time_keys=[k for k in self.__child.variables.keys() if k[:1] in ['I','D']]
		time_dim_keys=[(k[:1],k[1:6],k[6:]) for k in time_keys]
		split_names=[(k[:3],self.__groupsbyId[k[3:6]],self.__regionsbyId[k[6:]]) for k in spc_keys]
		named_keys=[self.__delim.join(k) for k in split_names]
		group_keys=list(set([k[0]+self.__delim+k[1]+self.__delim for k in split_names]))
		region_keys=list(set([k[0]+self.__delim*2+k[2] for k in split_names]))

		o3_keys=list(set([k[0]+self.__delim*2 for k in split_names]))
		self.variables=PseudoNetCDFVariables(self.__variables,spc_keys+time_dim_keys+named_keys+group_keys+region_keys+o3_keys)
	def __indices(self,keys):
		return [self.__child.__var_names__.index(k) for k in keys]

	def __var_id_names(self,var_name):
		if var_name in self.__child.variables.keys():
			keys=[var_name]
		elif var_name[:3] in ['O3N','O3V','NOX','VOC']:
			spc,group,region=var_name.split(self.__delim)
			espc=len(spc)
			egrp=3+len(group)
			ereg=6+len(region)
			keys=[k for k in self.__child.__var_names__ if spc=='' or spc==k[0:espc]]
			keys=[k for k in keys if group=='' or group==self.__groupsbyId[k[3:egrp]]]
			keys=[k for k in keys if region=='' or region==self.__regionsbyId[k[6:ereg]]]
		elif var_name[:1] in ['I','D']:
			spc,group,region=var_name.split(self.__delim)
			espc=1
			edayhr=1+5
			ereg=6+len(region)
			keys=[k for k in self.__child.__var_names__ if spc=='' or spc==k[0:espc]]
			keys=[k for k in keys if group=='' or group==self.__groupsbyId[k[1:edayhr]]]
			keys=[k for k in keys if region=='' or region==self.__regionsbyId[k[6:ereg]]]
		return keys

	def __var_nm_names(self,var_name):
		keys=self.__var_id_names(var_name)
		if var_name[:1] in ['I','D']:
			keys=[self.__delim.join([k[:1],self.__groupsbyId[k[1:6]],self.__regionsbyId[k[6:]]]) for k in keys]
		elif var_name[:3] in ['O3N','O3V','NOX','VOC']:
			keys=[self.__delim.join([k[:3],self.__groupsbyId[k[3:6]],self.__regionsbyId[k[6:]]]) for k in keys]
		
		return keys

	def __variables(self,key):
		var_id_names=self.__var_id_names(key)
		var_nm_names=self.__var_nm_names(key)
		
		var_indices=self.__indices(var_id_names)
		dimensions=('TSTEP','VAR','LAY','ROW','COL')
		outvals=self.__child.__memmap__[:,:,var_indices,:,:].swapaxes(1,2)
		v=PseudoNetCDFVariable(self,key,'f',dimensions,outvals)
		v.units='ppm'
		v.long_name=v.var_desc=key.ljust(16)
		v.VAR_NAMES=''.join([nm.ljust(16) for nm in var_id_names])
		v.VAR_NAME_DESCS=''.join([nm.ljust(16) for nm in var_nm_names])
		return v
