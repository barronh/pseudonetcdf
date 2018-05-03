__all__ = ['ncf2height_pressure', 'write_hgtprss']
__doc__ = """
.. _Write
:mod:`Write` -- CAMx height_pressure writer
============================================

.. module:: Write
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` writer for CAMx
              height_pressure files.
              for interface details
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""

# Site-Packages
from numpy import array

# This Package modules
from PseudoNetCDF.camxfiles.timetuple import timeadd, timerange
from PseudoNetCDF.camxfiles.FortranFileUtil import writeline
from PseudoNetCDF._getwriter import registerwriter


def ncf2height_pressure(ncffile, outpath, hght='HGHT', pres='PRES',
                        tflag='TFLAG'):
    outfile = open(outpath, 'wb')
    thp = zip(ncffile.variables[tflag][:, 0, :],
              ncffile.variables[hght],
              ncffile.variables[pres])
    for (d, t), h3d, p3d in thp:
        t = array(t.astype('>f') / 100, ndmin=1).astype('>f')
        d = array(d, ndmin=1).astype('>i')
        d = (d % (d // 100000 * 100000)).astype('>i')
        for i, (h2d, p2d) in enumerate(zip(h3d, p3d)):
            h2d = h2d.astype('>f')
            p2d = p2d.astype('>f')
            buf = array((h2d.size + 2) * 4, ndmin=1).astype('>i').tostring()
            outfile.write(buf + t.tostring() + d.tostring())
            h2d.tofile(outfile)
            outfile.write(buf)
            outfile.write(buf + t.tostring() + d.tostring())
            p2d.tofile(outfile)
            outfile.write(buf)
    outfile.flush()
    return outfile


registerwriter('camxfiles.height_pressure', ncf2height_pressure)
registerwriter('height_pressure', ncf2height_pressure)


def write_hgtprss(sdate, stime, time_step, vals):
    """
    Takes an iterable and some defining information
    and creates a CAMx read wind file

    sdate - integer start date
    stime - float start time
    time_step - float increment between times

    vals - array axes time,uv,xy,z
    """
    hp_string = ""
    edate, etime = timeadd((sdate, stime), (0, vals.shape[0] * time_step))

    tind = (2400, 24)[int(time_step % 2)]
    tenum = enumerate(timerange((sdate, stime),
                                timeadd((edate, etime), (0, 0)),
                                time_step,
                                tind))
    for i, (d, t) in tenum:
        for k in range(vals.shape[-1]):
            for hp in range(2):
                hpvals = [t, d]
                hpvals.extend(vals[i, hp, ..., k])
                hp_string += writeline(hpvals, "fi" + "f" * vals.shape[-2])
    return hp_string
