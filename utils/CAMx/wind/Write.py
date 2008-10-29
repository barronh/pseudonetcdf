HeadURL="$HeadURL: http://dawes.sph.unc.edu:8080/uncaqmlsvn/pyPA/utils/trunk/CAMxWrite.py $"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy: svnbarronh $"
__version__ = RevisionNum

#Distribution packages
from types import GeneratorType
import unittest
import struct,sys,os,operator
from warnings import warn
from tempfile import TemporaryFile as tempfile
import os,sys

#Site-Packages
from numpy import zeros,array,where,memmap,newaxis,dtype
try:
    from Scientific.IO.NetCDF import NetCDFFile as ncf
except:
    from pynetcdf import NetCDFFile as ncf

#This Package modules
from pyPA.utils.timetuple import timediff,timeadd,timerange
from pyPA.utils.util import cartesian,sliceit
from pyPA.utils.FortranFileUtil import OpenRecordFile,read_into,writeline,Int2Asc,Asc2Int
from pyPA.utils.sci_var import PseudoNetCDFFile, PseudoNetCDFVariable

def write_wind(sdate,stime,time_step,vals,lstagger=None):
    """
    Takes an iterable and some defining information
    and creates a CAMx read wind file
    
    sdate - integer start date
    stime - float start time
    time_step - float increment between times
    
    vals - array axes time,uv,xy,z
    """
    wind_string=""    
    edate,etime=timeadd((sdate,stime),(0,vals.shape[0]*time_step))
    
    for i,(d,t) in enumerate(timerange((sdate,stime),(edate,etime),time_step)):
        if lstagger!=None:
            wind_string+=writeline((t,d,lstagger),"fii")
        else:
            wind_string+=writeline((t,d),"fi")
            
        for k in range(vals.shape[-1]):
            for uv in range(2):
                wind_string+=writeline(vals[i,uv,...,k],"f"*vals.shape[-2])
        wind_string+=writeline((0,),"i")
    return wind_string

