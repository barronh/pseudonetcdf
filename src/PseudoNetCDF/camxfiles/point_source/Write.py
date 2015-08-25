__all__ = ['write_point']
__doc__ = """
.. _Write
:mod:`Write` -- CAMx point_source writer
============================================

.. module:: Write
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` writer for CAMx
              point_source files.  See PseudoNetCDF.sci_var.PseudoNetCDFFile 
              for interface details
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""

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
import numpy as np

#This Package modules
from PseudoNetCDF.camxfiles.timetuple import timediff,timeadd,timerange
from PseudoNetCDF.camxfiles.FortranFileUtil import OpenRecordFile,read_into,writeline,Int2Asc,Asc2Int
from PseudoNetCDF.sci_var import PseudoNetCDFFile, PseudoNetCDFVariable

_emiss_hdr_fmt=np.dtype(dict(names=['SPAD','name','note','itzon','nspec','ibdate','btime','iedate','etime','EPAD'],formats=['>i','(10,4)>S1','(60,4)>S1','>i','>i','>i','>f','>i','>f','>i']))

_grid_hdr_fmt=np.dtype(dict(names=['SPAD','plon','plat','iutm','xorg','yorg','delx','dely','nx','ny','nz','iproj','istag','tlat1','tlat2','rdum5','EPAD'],formats=['>i','>f','>f','>i','>f','>f','>f','>f','>i','>i','>i','>i','>i','>f','>f','>f','>i']))

_cell_hdr_fmt=np.dtype(dict(names=['SPAD','ione1','ione2','nx','ny','EPAD'],formats=['>i','>i','>i','>i','>i','>i']))

_time_hdr_fmt=np.dtype(dict(names=['SPAD','ibdate','btime','iedate','etime','EPAD'],formats=['>i','>i','>f','>i','>f','>i']))

    
_nstk_hdr_fmt=np.dtype(dict(names=['SPAD','ione','nstk','EPAD'],formats=['>i','>i','>i','>i']))
    
    
def ncf2point_source(ncffile, outpath):
    outfile = open(outpath, 'wb')
    emiss_hdr = np.zeros((1,), dtype = _emiss_hdr_fmt)
    emiss_hdr['SPAD'] = emiss_hdr['EPAD'] = emiss_hdr.dtype.itemsize - 8
    emiss_hdr[0]['name'][:, :] = ' '
    emiss_hdr[0]['name'][:, 0] = np.array([c for c in ncffile.NAME], dtype = '>S1')
    emiss_hdr[0]['note'][:, :] = ' '
    emiss_hdr[0]['note'][:, 0] = np.array([c for c in ncffile.NOTE], dtype = '>S1')
    gdtype = getattr(ncffile, 'GDTYPE', -999)
    emiss_hdr['itzon'][0] = ncffile.ITZON
    SPC_NAMES = [str(s) for s in np.array([c for c in ncffile.SPC_NAMES]).reshape(-1, 16)[:, :10].copy().view('>S10')[:,0]]
    nspec = len(SPC_NAMES)
    emiss_hdr['nspec'] = nspec

    grid_hdr = np.zeros((1,), dtype = _grid_hdr_fmt)
    grid_hdr['SPAD'] = grid_hdr['EPAD'] = grid_hdr.dtype.itemsize - 8
    grid_hdr['plon'] = ncffile.PLON
    grid_hdr['plat'] = ncffile.PLAT
    grid_hdr['iutm'][0] = ncffile.IUTM
    grid_hdr['xorg'] = ncffile.XORIG
    grid_hdr['yorg'] = ncffile.YORIG
    grid_hdr['delx'] = ncffile.XCELL
    grid_hdr['dely'] = ncffile.YCELL
    grid_hdr['nx'] = ncffile.NCOLS
    grid_hdr['ny'] = ncffile.NROWS
    grid_hdr['nz'] = ncffile.NLAYS
    grid_hdr['iproj'] = ncffile.CPROJ
    grid_hdr['tlat1'] = ncffile.TLAT1
    grid_hdr['tlat2'] = ncffile.TLAT2
    grid_hdr['istag'] = ncffile.ISTAG
    grid_hdr['rdum5'] = 100.
    
    cell_hdr = np.zeros((1,), dtype = _cell_hdr_fmt)
    cell_hdr['SPAD'] = cell_hdr['EPAD'] = cell_hdr.dtype.itemsize - 8
    cell_hdr['ione1'] = 1
    cell_hdr['ione2'] = 1
    cell_hdr['nx'] = ncffile.NCOLS
    cell_hdr['ny'] = ncffile.NROWS
    time_hdr = np.zeros(shape = (len(ncffile.dimensions['TSTEP']),), dtype = _time_hdr_fmt)
    time_hdr['SPAD'] = time_hdr['EPAD'] = time_hdr.dtype.itemsize - 8
    date, time = ncffile.variables['TFLAG'][:, 0].T
    edate, etime = ncffile.variables['ETFLAG'][:, 0].T
    time = time.astype('>f') / 10000.
    date = date%(date/100000*100000)
    etime = etime.astype('>f') / 10000.
    edate = edate%(edate/100000*100000)
    time_hdr['ibdate'] = date
    time_hdr['btime'] = time
    time_hdr['iedate'] = edate
    time_hdr['etime'] = etime
    #time_hdr['iedate'] += time_hdr['etime'] // 24
    #time_hdr['etime'] -= (time_hdr['etime'] // 24) * 24
    
    
    emiss_hdr['ibdate'] = time_hdr[0]['ibdate']
    emiss_hdr['btime'] = time_hdr[0]['btime']
    emiss_hdr['iedate'] = ncffile.EDATE
    emiss_hdr['etime'] = ncffile.ETIME
    
    nspec = emiss_hdr['nspec'][0]
    _spc_hdr_fmt=np.dtype(dict(names=['SPAD','NAME','EPAD'], formats = ['>i', '(%d,10,4)>S1' % nspec, '>i']))
    spc_hdr = np.zeros((1,), dtype = _spc_hdr_fmt)
    spc_hdr['NAME'][:] = ' '
    spc_hdr['SPAD'] = spc_hdr['EPAD'] = spc_hdr.dtype.itemsize - 8
    for si, spc in enumerate(SPC_NAMES):
        spc_hdr[0]['NAME'][si,:,0] = np.array([c for c in spc.ljust(10)], dtype = '>S1')
    
    nstk_hdr = np.ones((1,), dtype = _nstk_hdr_fmt)
    nstk_hdr['SPAD'] = nstk_hdr['EPAD'] = nstk_hdr.dtype.itemsize - 8
    nstk_hdr['nstk'] = len(ncffile.dimensions['NSTK'])
    
    nstk=nstk_hdr['nstk'][0]
    #(xstk(n),ystk(n),hstk(n),dstk(n),tstk(n),vstk(n),n=1,nstk)
    stk_prop_fmt = np.dtype(dict(names = ['xstk', 'ystk', 'hstk', 'dstk', 'tstk', 'vstk'], formats = ['>f', '>f', '>f', '>f', '>f', '>f']))
    stk_prop = np.ones((1,), dtype = np.dtype([('SPAD', '>i', 1), ('DATA', stk_prop_fmt, nstk), ('EPAD', '>i', 1)]))
    stk_prop['SPAD'] = stk_prop['EPAD'] = stk_prop.dtype.itemsize - 8
    stk_prop['DATA']['xstk'] = ncffile.variables['XSTK'][:]
    stk_prop['DATA']['ystk'] = ncffile.variables['YSTK'][:]
    stk_prop['DATA']['hstk'] = ncffile.variables['HSTK'][:]
    stk_prop['DATA']['dstk'] = ncffile.variables['DSTK'][:]
    stk_prop['DATA']['tstk'] = ncffile.variables['TSTK'][:]
    stk_prop['DATA']['vstk'] = ncffile.variables['VSTK'][:]
    
    # Start writing
    emiss_hdr.tofile(outfile)
    grid_hdr.tofile(outfile)
    cell_hdr.tofile(outfile)
    spc_hdr.tofile(outfile)
    nstk_hdr.tofile(outfile)
    spc_names = [str(np.char.strip(spc[:, 0].copy().view('>S10')[0])) for spc in spc_hdr[0]['NAME']]
    stk_prop.tofile(outfile)
    for di, (d, t) in enumerate(ncffile.variables['TFLAG'][:, 0]):
        time_hdr[di].tofile(outfile)
        np.array([8, 1, nstk, 8], dtype = '>i').tofile(outfile)
        #(idum,idum,kcell(n),flow(n),plmht(n),n=1,nstk)
        time_prop_fmt = np.dtype(dict(names = ['ione1', 'ione2', 'kcell', 'flow', 'plmht'], formats = ['>i', '>i', '>i', '>f', '>f']))
        props = np.ones((1,), dtype = np.dtype([('SPAD', '>i', 1), ('DATA', time_prop_fmt, nstk), ('EPAD', '>i', 1)]))
        props['SPAD'] = props['EPAD'] = props.dtype.itemsize - 8
        props['DATA']['ione1'] = ncffile.variables['IONE'][di]
        props['DATA']['ione2'] = ncffile.variables['ITWO'][di]
        props['DATA']['kcell'] = ncffile.variables['KCELL'][di]
        props['DATA']['flow'] = ncffile.variables['FLOW'][di]
        props['DATA']['plmht'] = ncffile.variables['PLMHT'][di]
        props.tofile(outfile)
        for spc_key, spc_name in zip(spc_names, spc_hdr[0]['NAME']):
            var = ncffile.variables[spc_key]
            data = var[di].astype('>f')
            buf = np.array([4+40+data.size*4]).astype('>i')
            buf.tofile(outfile)
            np.array([1]).astype('>i').tofile(outfile)
            spc_name.tofile(outfile)
            data.tofile(outfile)
            buf.tofile(outfile)
    outfile.flush()
    return outfile

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