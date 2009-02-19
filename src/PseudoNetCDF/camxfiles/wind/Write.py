__all__ = ['write_wind']
__doc__ = """
.. _Write
:mod:`Write` -- CAMx wind writer
============================================

.. module:: Write
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` writer for CAMx
              wind files.  See PseudoNetCDF.sci_var.PseudoNetCDFFile 
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

#This Package modules
from PseudoNetCDF.camxfiles.timetuple import timeadd,timerange
from PseudoNetCDF.camxfiles.FortranFileUtil import writeline

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

