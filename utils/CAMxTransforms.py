try:
    from Scientific.IO.NetCDF import NetCDFFile as ncf
except:
    from pynetcdf import NetCDFFile as ncf
from sci_var import PseudoNetCDFFile, PseudoIOAPIVariable
from CAMx.wind.Transforms import wind_center_time_cell
from CAMx.height_pressure.Transforms import height_pressure_center_time_plus, \
                                            height_pressure_center_time, \
                                            height_pressure_plus
from CAMx.temperature.Transforms import temperature_center_time
from CAMx.vertical_diffusivity.Transforms import vertical_diffusivity_center_time
from CAMx.humidity.Transforms import humidity_center_time
from CAMx.cloud_rain.Transforms import cloud_rain_center_time_plus, \
                                       cloud_rain_center_time, \
                                       cloud_rain_plus
from CAMxFiles import *
from MetaNetCDF import *
from CAMxFiles import point_source as reg_point_source, \
                      gridded_emissions as reg_gridded_emissions
from sci_var import PseudoNetCDFFile, \
                    PseudoNetCDFVariables, \
                    PseudoNetCDFVariable


__all__ = ['point_source_newvarnames', 'pypass_camx_met_master', 'camx_pa_master', 'hght2dz', 'hght2zh', 'pypass_camx_emiss_master']

#==================================================================

#                                                             point_source 

class point_source_newvarnames(PseudoNetCDFFile):
    __childclass__=reg_point_source
    __newvars__=['ptCO','ptNO','ptNO2','ptALD2','ptETH','ptFORM','ptISOP','ptNR','ptOLE','ptPAR','ptTOL','ptXYL', 'ptNH3', 'ptSO2', 'ptSULF', 'ptPEC','ptPNO3','ptPOA','ptPSO4','ptETOH','ptMEOH']
    __oldvars__=['XSTK','YSTK','HSTK']  
    def __init__(self,rffile):
        self.__child=self.__childclass__(rffile)
        self.dimensions=self.__child.dimensions
        self.variables=PseudoNetCDFVariables(self.__variables__,self.__newvars__ + self.__oldvars__)
    def __variables__(self,key):
        if key in self.__newvars__:
           val = self.__child.variables[key[2:]]
           var=PseudoNetCDFVariable(self,key,'f',('TSTEP','STK'),values=val)
        elif key in self.__oldvars__:
           val = self.__child.variables[key]
           var=PseudoNetCDFVariable(self,key,'f',('TSTEP','STK'),values=val)
#          return getattr(self,var)
        return var
                
#==================================================================

def pypass_camx_met_master(wind_path,hp_path,temp_path,kv_path,hum_path,cr_path,rows,cols,endhour=True):
     windf=wind_center_time_cell(wind_path,rows,cols,outunit='km/h',endhour=endhour,forcestaggered=False)
     hpf=height_pressure_center_time_plus(hp_path,rows,cols,outunit={'HGHT':'km', 'PRES':'hPA'},endhour=endhour)
     tempf=temperature_center_time(temp_path,rows,cols,outunit={'AIRTEMP':'deg_F','SURFTEMP':'deg_F'},endhour=endhour)
     kvf=vertical_diffusivity_center_time(kv_path,rows,cols,outunit={'KV':'m**2/s'},endhour=endhour)
     humf=humidity_center_time(hum_path,rows,cols,outunit={'HUM':'ppm'},endhour=endhour)
     crf=cloud_rain_center_time_plus(cr_path,rows,cols,outunit={'CLOUD':'g/m**3','RAIN':'g/m**3','SNOW':'g/m**3','GRAUPEL':'g/m**3','PRECIP':'g/m**3','PRECIPRATE':'mm/h','COD':'None'},endhour=endhour)
     return file_master([windf,hpf,tempf,kvf,humf,crf])

def camx_pa_master(paths_and_readers,tslice=slice(None),kslice=slice(None),jslice=slice(None),islice=slice(None)):
    """
    CAMx PA Master presents a single interface for CAMx IPR, IRR, and shape definitions.
    If the shape is not defined, CAMx PA Master can provide a default shape horiztonally
    bound by the domain and vertically bound by the planetary boundary layer.  The
    planetary boundary layer is diagnosed by from the vertical diffusivity and vertical
    layer structure.
    
    paths_and_readers - iterable of iterables (n x 2) where each element of the
                        primary iterable is an iterable containing a file path 
                        and a reader for that path.  The reader is expected to
                        present the Scientific.IO.NetCDF.NetCDFFile file inter-
                        face.
    
    optional:
       tslice - slice object used to window the time period from the 
                instantaneous meteorological files to match the PA domain
       kslice - same as tslice, but for layers
       jslice - same as tslice, but for rows
       islice -  - same as tslice, but for columns
              
    """
    def defaultshape(self):
        from pyPA.pappt.kvextract import tops2shape, vertcamx
        
        old_shape=[i for i in self.variables['UCNV_O3'].shape]
        new_shape=[i for i in self.variables['UCNV_O3'].shape]
        new_shape[0]+=1
        new_shape=zeros(tuple(new_shape),dtype='bool')
        new_shape[1:,:,:,:]=tops2shape(vertcamx(CenterTime(self.variables['KV']),CenterTime(self.variables['HGHT']))[tslice,jslice,islice],old_shape)
        new_shape[0,:,:,:]=new_shape[1,:,:,:]
        return PseudoIOAPIVariable(self,'DEFAULT_SHAPE','i',('TSTEP', 'LAY', 'ROW', 'COL'),values=new_shape,units='on/off')

    # Create a list of opened files
    files=[eval(r)(p) for p,r in paths_and_readers]
    
    # Create a file master object from files
    master_file=file_master(files)
    
    ## Add derived variables
    
    # Rename AVOL_O3 to VOL -- simplicity
    #   The choice of O3 was arbitrary.  All AVOL 
    #   values (e.g. AVOL_O3, AVOL_NO2, etc) are equal.
    master_file.addMetaVariable('VOL',lambda self: PseudoIOAPIVariable(self,'VOL','f',('TSTEP','LAY','ROW','COL'),values=self.variables['AVOL_O3'],units='m**3'))
    
    # Calculate AIRMOLS from IPR outputs
    #   UCNV [=] m**3/mol_{air}
    #   AVOL [=] m**3
    #   AIRMOLS = AVOL/UCNV [=] mol_air
    master_file.addMetaVariable('AIRMOLS',lambda self: PseudoIOAPIVariable(self,'AIRMOLS','f',('TSTEP','LAY','ROW','COL'),values=self.variables['AVOL_O3'][:]/self.variables['UCNV_O3'][:],units='moles**-1'))

    # Calculate INVAIRMOLS from IPR outputs
    #   UCNV [=] m**3/mol_{air}
    #   AVOL [=] m**3
    #   AIRMOLS = UCNV/AVOL [=] 1/mol_air
    master_file.addMetaVariable('INVAIRMOLS',lambda self: PseudoIOAPIVariable(self,'INVAIRMOLS','f',('TSTEP','LAY','ROW','COL'),values=self.variables['UCNV_O3']/self.variables['AVOL_O3'],units='moles**-1'))

    # Calculate the well mixed region of the troposphere based
    # on vertical diffusivity and layer structure.
    master_file.addMetaVariable('DEFAULT_SHAPE',lambda self: defaultshape(self))
    
    return master_file


def hght2dz(layer_tops):
    dz=layer_tops.copy()
    dz[:,1:,:,:]-=layer_tops[:,:-1,:,:]
    return dz

def hght2zh(layer_tops):
    return layer_tops-.5*dz

def pypass_camx_emiss_master(pointemiss_path,gridemiss_path,rows,cols,endhour=True):
    pointemiss=point_source_newvarnames(pointemiss_path)
    griddedemiss=reg_gridded_emissions(gridemiss_path)
    return file_master([pointemiss,griddedemiss])
    