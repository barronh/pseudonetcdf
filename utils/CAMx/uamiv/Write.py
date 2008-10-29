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

def ncf2uamiv(ncffile,outpath):
	outfile=file(outpath,'wb')
	for val3dw in zip(ncffile.variables.values()):
		t=array(t.astype('>f')/10000,ndmin=1).astype('>f')
		d=array(d,ndmin=1).astype('>i')
		d=(d%(d/100000*100000)).astype('>i')
		for i,(h2d,p2d) in enumerate(zip(h3d,p3d)):
			h2d=h2d.astype('>f')
			p2d=p2d.astype('>f')
			buf=array((h2d.size+2)*4,ndmin=1).astype('>i').tostring()
			outfile.write(buf+t.tostring()+d.tostring()+h2d.tostring()+buf)
			outfile.write(buf+t.tostring()+d.tostring()+p2d.tostring()+buf)

def write_emissions_ncf(infile,outfile):
    from operator import concat
    #initialize hdr_fmts with species count
    hdr_fmts=["10i60i3ifif","ffiffffiiiiifff","iiii","10i"*len(infile.variables.keys())]
    hdrlines=[]

    hdrlines.append(reduce(concat,[Asc2Int(s) for s in [infile.name, infile.note]])+[infile.ione, len(infile.variables.keys()),infile.start_date,infile.start_time,infile.end_date,infile.end_time])

    hdrlines.append([infile.rdum, infile.rdum, infile.iutm, infile.xorg, infile.yorg, infile.delx, infile.dely, infile.dimensions['COL'], infile.dimensions['ROW'], infile.dimensions['LAY'], infile.idum, infile.idum, infile.rdum, infile.rdum, infile.rdum])
    
    hdrlines.append([infile.ione,infile.ione,infile.dimensions['COL'],infile.dimensions['ROW']])
    hdrlines.append(reduce(concat,[Asc2Int(s.ljust(10)) for s in infile.variables.keys()]))

    for d,h in zip(hdrlines,hdr_fmts):
        outfile.write(writeline(d,h))
        
    for ti,(d,t) in enumerate(infile.timerange()):
        ed,et=timeadd((d,t),(0,infile.time_step))
        outfile.write(writeline((d,t,ed,et),'ifif'))
        for spc in infile.variables.keys():
            for k in range(infile.dimensions['LAY']):
                outfile.write(writeline([infile.ione]+Asc2Int(spc.ljust(10))+infile.variables[spc][ti,:,:,k].transpose().ravel().tolist(),'11i'+infile.cell_count*'f'))

def write_emissions(start_date,start_time,time_step,hdr,vals):
    #initialize hdr_fmts with species count
    hdr_fmts=["10i60i3ifif","ffiffffiiiiifff","iiii","10i"*len(hdr[-1])]
    
    #initialize output variable
    emis_string=''
    
    #Internalize header
    hdr=[h for h in hdr]
    species=hdr[-1]
    
    #Check start_date and start_time
    if tuple(hdr[0][4:6])!=(start_date,start_time):
        print >>sys.stderr,"Header doesn't match start date/time"
       
    #Change name and note
    hdr[0]=list(hdr[0][0])+list(hdr[0][1])+list(hdr[0][2:])
    
    #Change species names to array of characters
    hdr[-1]=reduce(operator.concat,[Asc2Int(s) for s in hdr[-1]])
    
    #for each item in the header, write it to output
    for h,f in zip(hdr,hdr_fmts):
        emis_string+=writeline(h,f)

    #create value format
    cells=vals.shape[2]*vals.shape[3]
    valfmt='i10i'+('f'*cells)
    
    #Get end date
    (end_date,end_time)=timeadd((start_date,start_time),(0,time_step*vals.shape[0]))
    
    if tuple(hdr[0][6:])!=(end_date,end_time):
        print >>sys.stderr,"Header doesn't match end date/time"
    
    #Write out values
    for ti,(d,t) in enumerate(timerange((start_date,start_time),(end_date,end_time),time_step)):
        ed,et=timeadd((d,t),(0,time_step))
        emis_string+=writeline((d,t,ed,et),'ifif')
        for si,spc in enumerate(species):
            for k in range(vals.shape[-1]):
                #Dummy variable,spc characters and values flattened
                temp=[1]
                temp.extend(Asc2Int(spc))
                spcvals=vals[ti,si,...,...,k]
                spcvals=spcvals.transpose().ravel()
                temp.extend(spcvals)
                emis_string+=writeline(temp,valfmt)
    
    return emis_string        

