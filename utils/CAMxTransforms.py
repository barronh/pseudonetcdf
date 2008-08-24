from MetaNetCDF import file_master, \
                       time_avg_new_unit
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
from CAMxFiles import point_source as reg_point_source, \
                      gridded_emissions as reg_gridded_emissions
from sci_var import PseudoNetCDFFile, \
                    PseudoNetCDFVariables, \
                    PseudoNetCDFVariable
                                 
#==================================================================

#                                                              point_source 

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
           var=PseudoNetCDFVariable(self,key,'f',('TSTEP','STK'),val)
        elif key in self.__oldvars__:
           val = self.__child.variables[key]
           var=PseudoNetCDFVariable(self,key,'f',('TSTEP','STK'),val)
#           return getattr(self,var)
        return var
				
#==================================================================

def pypass_camx_met_master(wind_path,hp_path,temp_path,kv_path,hum_path,cr_path,rows,cols,endhour=True):
     windf=wind_center_time_cell(wind_path,rows,cols,outunit='km/h',endhour=endhour,forcestaggered=False)
     hpf=height_pressure_center_time_plus(hp_path,rows,cols,outunit={'HGHT':'km', 'PRES':'hPA'},endhour=endhour)
     tempf=temperature_center_time(temp_path,rows,cols,outunit={'AIRTEMP':'deg_F','SURFTEMP':'deg_F'},endhour=endhour)
     kvf=vertical_diffusivity_center_time(kv_path,rows,cols,outunit={'KV':'m**2/s'},endhour=endhour)
     humf=humidity_center_time(hum_path,rows,cols,outunit={'HUM':'ppm'},endhour=endhour)
     crf=cloud_rain_center_time_plus(cr_path,rows,cols,outunit={'CLOUD':'g/m**3','RAIN':'g/m**3','SNOW':'g/m**3','GRAUPEL':'g/m**3','PRECIP':'g/m**3','PRECIP_RATE':'mm/h','COD':'None'},endhour=endhour)
     return file_master([windf,hpf,tempf,kvf,humf,crf])

def pypass_camx_emiss_master(pointemiss_path,gridemiss_path,rows,cols,endhour=True):
    pointemiss=point_source_newvarnames(pointemiss_path)
    griddedemiss=reg_gridded_emissions(gridemiss_path)
    return file_master([pointemiss,griddedemiss])