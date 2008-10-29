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

def write_point(start_date,start_time,time_step,hdr,vals):
    #Internalize header
    hdr=[h for h in hdr]
    species=hdr[-4]
    nstk=hdr[-3][1]
    timeprops=hdr.pop()
    
    #initialize hdr_fmts with species count
    hdr_fmts=["10i60i3ifif","ffiffffiiiiifff","iiii","10i"*len(species),"ii","ffffff"*nstk]
    
    #initialize output variable
    pt_string=''
    
    #Change name and note
    hdr[0]=list(hdr[0][0])+list(hdr[0][1])+list(hdr[0][2:])
    
    #Reducing stk props
    stkprops=hdr[-1]
    hdr[-1]=[]
    for stk in stkprops:
        hdr[-1].extend(stk)
    stk_time_prop_fmt="iiiff"*nstk
    
    stk_time_props=[]
    for time in timeprops:
        stk_time_props.append([])
        for stk in time:
            stk_time_props[-1].extend(stk)
            
    #Change species names to array of characters
    hdr[-3]=reduce(operator.concat,[Asc2Int(s) for s in hdr[-3]])
    
    #for each item in the header, write it to output
    for i,(h,f) in enumerate(zip(hdr,hdr_fmts)):
        pt_string+=writeline(h,f,False)

    #create value format
    valfmt='i10i'+('f'*nstk)
    
    #Get end date
    (end_date,end_time)=timeadd((start_date,start_time),(0,time_step*vals.shape[0]))
    
    #Write out values
    for ti,(d,t) in enumerate(timerange((start_date,start_time),(end_date,end_time),time_step)):
        ed,et=timeadd((d,t),(0,time_step))
        pt_string+=writeline((d,t,ed,et),'ifif',False)
        pt_string+=writeline((1,nstk),'ii',False)
        pt_string+=writeline(stk_time_props[ti],stk_time_prop_fmt,False)
        for si,spc in enumerate(species):
            #Dummy variable,spc characters and values flattened
            temp=[1]
            temp.extend(Asc2Int(spc))
            temp.extend(vals[ti,si,...])
            pt_string+=writeline(temp,valfmt,False)
    
    return pt_string