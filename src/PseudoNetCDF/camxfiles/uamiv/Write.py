__all__ = ['ncf2uamiv', 'write_emissions_ncf', 'write_emissions']
__doc__ = """
.. _Write
:mod:`Write` -- CAMx uamiv  writer
============================================

.. module:: Write
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` writer for CAMx
              uamiv files.  See PseudoNetCDF.sci_var.PseudoNetCDFFile 
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
from numpy import zeros,array,where,memmap,newaxis,dtype
import numpy as np

#This Package modules
from PseudoNetCDF.camxfiles.timetuple import timediff,timeadd,timerange
from PseudoNetCDF.camxfiles.FortranFileUtil import OpenRecordFile,read_into,writeline,Int2Asc,Asc2Int
from PseudoNetCDF.sci_var import PseudoNetCDFFile, PseudoNetCDFVariable

_emiss_hdr_fmt=np.dtype(dict(names=['SPAD','name','note','itzon','nspec','ibdate','btime','iedate','etime','EPAD'],formats=['>i','(10,4)>S1','(60,4)>S1','>i','>i','>i','>f','>i','>f','>i']))

_grid_hdr_fmt=np.dtype(dict(names=['SPAD','plon','plat','iutm','xorg','yorg','delx','dely','nx','ny','nz','iproj','istag','tlat1','tlat2','rdum5','EPAD'],formats=['>i','>f','>f','>i','>f','>f','>f','>f','>i','>i','>i','>i','>i','>f','>f','>f','>i']))

_cell_hdr_fmt=np.dtype(dict(names=['SPAD','ione1','ione2','nx','ny','EPAD'],formats=['>i','>i','>i','>i','>i','>i']))

_time_hdr_fmt=np.dtype(dict(names=['SPAD','ibdate','btime','iedate','etime','EPAD'],formats=['>i','>i','>f','>i','>f','>i']))

_spc_fmt=np.dtype("(10,4)>S1")

def ncf2uamiv(ncffile, outpath):
    """
    ncf2uamiv converts a ncffile to a uamiv file
    
    Parameters
    ----------
    ncffile : PseudoNetCDF-like object
            Must be open and have all properties of a uamiv
            file as produced by PseudoNetCDF.camxfiles
            1) IOAPI conventions
            2) UAMIV specific properties: PLON PLAT IUTM CPROJ 
                                            TLAT1 TLAT2 ISTAG
    outpath : string
            path to create a uamiv file output
    
    Returns
    -------
    outfile : file object
            File object is in write binary mode and has the
            uamiv file as its contents

    Examples
    --------
        $ # test.uamiv must be a file in uamiv format
        $ pncgen -f uamiv test.uamiv test.uamiv.nc
        $ python -c "
        ncpath = 'test.uamiv.nc'
        uamivpath = 'test.uamiv.nc.uamiv'
        from netCDF4 import Dataset
        from PseudoNetCDF.camxfiles.uamiv.Write import ncf2uamiv
        inncf = Dataset(ncpath)
        ncf2uamiv(inncf, uamivpath)
        "
        # If uamivpath does not include EOD, the diff will be perfect
        diff test.uamiv test.uamiv.nc.uamiv

        # If uamivpath includes EOD, the diff may yield difference.
        # The ONLY difference will be time flags for the end of a day.
        # Most time software agrees that there is no such thing as 2002154T24:00.
        # Some CAMx files, however, have a 2400 time. 
        # PseudoNetCDF interprets this equivalently as 2002155T00:00
        $ python -c "from PseudoNetCDF.camxfiles.Memmaps import uamiv
        import numpy as np
        old = uamiv('test.uamiv')
        new = uamiv('test.uamiv.nc.uamiv')
        for k in old.variables.keys():
            check = (old.variables[k][...] == new.variables[k][...])
            if not check.all():
                print k
                if k == 'ETFLAG':
                    diffidx = np.where(~check)[:2]
                else:
                    diffidx = np.where(~check)
                print old.variables[k][diffidx]
                print new.variables[k][diffidx]

    """
    
    emiss_hdr = np.zeros(shape = (1,), dtype = _emiss_hdr_fmt)
    emiss_hdr[0]['name'][:, :] = ' '
    emiss_hdr[0]['name'][:, 0] = np.array(ncffile.NAME, dtype = '>c')
    emiss_hdr[0]['note'][:, :] = ' '
    emiss_hdr[0]['note'][:, 0] = np.array(ncffile.NOTE, dtype = '>c')
    gdtype = getattr(ncffile, 'GDTYP', -999)
    emiss_hdr['itzon'][0] = ncffile.ITZON
    nspec = len(ncffile.dimensions['VAR'])
    emiss_hdr['nspec'] = nspec
    
    NCOLS = len(ncffile.dimensions['COL'])
    NROWS = len(ncffile.dimensions['ROW'])
    NLAYS = len(ncffile.dimensions['LAY'])
    grid_hdr = np.zeros(shape = (1,), dtype = _grid_hdr_fmt)
    grid_hdr['SPAD'] = grid_hdr.itemsize - 8
    grid_hdr['plon'] = ncffile.PLON
    grid_hdr['plat'] = ncffile.PLAT
    grid_hdr['iutm'][0] = ncffile.IUTM
    grid_hdr['xorg'] = ncffile.XORIG
    grid_hdr['yorg'] = ncffile.YORIG
    grid_hdr['delx'] = ncffile.XCELL
    grid_hdr['dely'] = ncffile.YCELL
    grid_hdr['nx'] = NCOLS
    grid_hdr['ny'] = NROWS
    grid_hdr['nz'] = NLAYS
    grid_hdr['iproj'] = ncffile.CPROJ
    grid_hdr['tlat1'] = ncffile.TLAT1
    grid_hdr['tlat2'] = ncffile.TLAT2
    grid_hdr['istag'] = ncffile.ISTAG
    grid_hdr['rdum5'] = 0.
    grid_hdr['EPAD'] = grid_hdr.itemsize - 8

    cell_hdr = np.zeros(shape = (1,), dtype = _cell_hdr_fmt)
    cell_hdr['SPAD'] = cell_hdr.itemsize - 8
    cell_hdr['ione1'] = 1
    cell_hdr['ione2'] = 1
    cell_hdr['nx'] = NCOLS
    cell_hdr['ny'] = NROWS
    cell_hdr['EPAD'] = cell_hdr.itemsize - 8

    time_hdr = np.zeros(shape = (len(ncffile.dimensions['TSTEP']),), dtype = _time_hdr_fmt)
    time_hdr['SPAD'] = 16
    time_hdr['EPAD'] = 16
    date, time = ncffile.variables['TFLAG'][:, 0].T
    time = time.astype('>f') / 10000.
    date = date%(date/100000*100000)
    time_hdr['ibdate'] = date
    time_hdr['btime'] = time
    time_hdr['iedate'] = date
    time_hdr['etime'] = time + 1.
    time_hdr['iedate'] += time_hdr['etime'] // 24
    time_hdr['etime'] -= (time_hdr['etime'] // 24) * 24
    emiss_hdr['ibdate'] = time_hdr['ibdate'][0]
    emiss_hdr['btime'] = time_hdr['btime'][0]
    emiss_hdr['iedate'] = time_hdr['iedate'][0]
    emiss_hdr['etime'] = time_hdr['etime'][0]
    emiss_hdr['SPAD'] = _emiss_hdr_fmt.itemsize - 8
    emiss_hdr['EPAD'] = _emiss_hdr_fmt.itemsize - 8
    
    spc_hdr = np.zeros(shape = (1,), dtype = dict(names = ['SPAD1', 'DATA', 'EPAD1'], formats = ['>i', np.dtype("(%d,10,4)>S1" % nspec), '>i']))
    spc_hdr['SPAD1'] = nspec * 40
    spc_hdr['EPAD1'] = nspec * 40
    spc_names = np.array(getattr(ncffile, 'VAR-LIST'), dtype = '>c').reshape(-1, 16)[:, :10].copy()
    spc_hdr[0]['DATA'][:] = ' '
    spc_hdr[0]['DATA'][:, :, 0] = spc_names
    spc_names = spc_names.view('>S10')
    nz = len(ncffile.dimensions['LAY'])
    outfile = open(outpath, 'wb')
    emiss_hdr.tofile(outfile)
    grid_hdr.tofile(outfile)
    cell_hdr.tofile(outfile)
    spc_hdr.tofile(outfile)
    for di, (d, t) in enumerate(ncffile.variables['TFLAG'][:, 0]):
        time_hdr[di].tofile(outfile)
        for spc_key, spc_name in zip(spc_names, spc_hdr[0]['DATA']):
            for zi in range(nz):
                var = ncffile.variables[str(np.char.strip(spc_key[0]))]
                data = var[di, zi].astype('>f')
                buf = np.array(4+40+data.size*4).astype('>i')
                buf.tofile(outfile)
                np.array(1).astype('>i').tofile(outfile)
                spc_name.tofile(outfile)
                data.tofile(outfile)
                buf.tofile(outfile)
    outfile.flush()
    return outfile

def write_emissions_ncf(infile,outfile):
    from operator import concat
    #initialize hdr_fmts with species count
    hdr_fmts=["10i60i3ifif","ffiffffiiiiifff","iiii","10i"*len(infile.variables.keys())]
    hdrlines=[]

    hdrlines.append(reduce(concat,[Asc2Int(s) for s in [infile.name, infile.note]])+[infile.ione, len(infile.variables.keys()),infile.start_date,infile.start_time,infile.end_date,infile.end_time])

    hdrlines.append([infile.rdum, infile.rdum, infile.iutm, infile.xorg, infile.yorg, infile.delx, infile.dely, len(infile.dimensions['COL']), len(infile.dimensions['ROW']), len(infile.dimensions['LAY']), infile.idum, infile.idum, infile.rdum, infile.rdum, infile.rdum])
    
    hdrlines.append([infile.ione,infile.ione,len(infile.dimensions['COL']),len(infile.dimensions['ROW'])])
    hdrlines.append(reduce(concat,[Asc2Int(s.ljust(10)) for s in infile.variables.keys()]))

    for d,h in zip(hdrlines,hdr_fmts):
        outfile.write(writeline(d,h))
        
    for ti,(d,t) in enumerate(infile.timerange()):
        ed,et=timeadd((d,t),(0,infile.time_step))
        outfile.write(writeline((d,t,ed,et),'ifif'))
        for spc in infile.variables.keys():
            for k in range(len(infile.dimensions['LAY'])):
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

