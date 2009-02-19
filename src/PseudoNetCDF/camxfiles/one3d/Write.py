__all__ = ['ncf2one3d']
__doc__ = """
.. _Write
:mod:`Write` -- CAMx one3d writer
============================================

.. module:: Write
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` writer for CAMx generic
              1 3D variable files.  See PseudoNetCDF.sci_var.PseudoNetCDFFile 
              for interface details
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""

from numpy import array

def ncf2one3d(ncffile,outpath,key='KV',tflag='TFLAG'):
    outfile=file(outpath,'wb')
    for (d,t),v3d in zip(ncffile.variables[tflag][:,0,:],ncffile.variables[key]):
        t=array(t.astype('>f')/10000,ndmin=1).astype('>f')
        d=array(d,ndmin=1).astype('>i')
        d=(d%(d/100000*100000)).astype('>i')
        for v2d in v3d:
            v2d=v2d.astype('>f')
            buf=array((v2d.size+2)*4,ndmin=1).astype('>i').tostring()
            outfile.write(buf+t.tostring()+d.tostring()+v2d.tostring()+buf)