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

import numpy as np

# This Package modules
from PseudoNetCDF.camxfiles.timetuple import timeadd, timerange
from PseudoNetCDF.camxfiles.FortranFileUtil import writeline
from PseudoNetCDF._getwriter import registerwriter

"""
hour,idate,lstagger
Loop from 1 to nlay layers:
    ((uwind(i,j,k),i=1,nx),j=1,ny)
    ((vwind(i,j,k),i=1,nx),j=1,ny)
((dummy,i=1,nx),j=1,ny)
"""


def ncf2wind(ncffile, outpath, tflag='TFLAG'):
    outfile = open(outpath, 'wb')
    varkeys = [vk for vk in ['U', 'V'] if vk in ncffile.variables.keys()]
    nzcl = len(ncffile.dimensions['LAY'])
    for di, (d, t) in enumerate(ncffile.variables[tflag][:, 0, :]):
        t = np.array(t // 100, ndmin=1).astype('>f')
        d = np.array(d, ndmin=1).astype('>i')
        d = (d % (d // 100000 * 100000)).astype('>i')
        lstag = ncffile.LSTAGGER
        buf = np.array([12], dtype='>i').tobytes()
        outfile.write(buf + t.tobytes() + d.tobytes() +
                      lstag.tobytes() + buf)
        for zi in range(nzcl):
            for varkey in varkeys:
                vals = ncffile.variables[varkey][di, zi].astype('>f')
                buf = np.array([(vals.size) * 4],
                               ndmin=1).astype('>i').tobytes()
                outfile.write(buf)
                vals.tofile(outfile)
                outfile.write(buf)
        vals = np.array(0, dtype='>i')
        buf = np.array([(vals.size) * 4], ndmin=1).astype('>i').tobytes()
        outfile.write(buf)
        vals.astype('>f').tofile(outfile)
        outfile.write(buf)
    outfile.flush()
    return outfile


registerwriter('camxfiles.wind', ncf2wind)
registerwriter('wind', ncf2wind)


def write_wind(sdate, stime, time_step, vals, lstagger=None):
    """
    Takes an iterable and some defining information
    and creates a CAMx read wind file

    sdate - integer start date
    stime - float start time
    time_step - float increment between times

    vals - array axes time,uv,xy,z
    """
    wind_string = ""
    edate, etime = timeadd((sdate, stime), (0, vals.shape[0] * time_step))

    for i, (d, t) in enumerate(timerange((sdate, stime), (edate, etime),
                                         time_step)):
        if lstagger is not None:
            wind_string += writeline((t, d, lstagger), "fii")
        else:
            wind_string += writeline((t, d), "fi")

        for k in range(vals.shape[-1]):
            for uv in range(2):
                wind_string += writeline(vals[i, uv, ..., k],
                                         "f" * vals.shape[-2])
        wind_string += writeline((0,), "i")
    return wind_string
