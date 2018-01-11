from netCDF4 import Dataset
from ..core._files import PseudoNetCDFFile
class ioapi_base(PseudoNetCDFFile):
    def __getattributte__(self, *args, **kwds):
        return getattr(self, *args, **kwds)
    @classmethod
    def isMine(self, path):
        return False
    
    def _newlike(self):
        if isinstance(self, PseudoNetCDFFile):
            outt = ioapi_base
            outf = outt.__new__(outt)
        else:
            outf = PseudoNetCDFFile()
        return outf 
    def setncatts(self, attdict):
        PseudoNetCDFFile.setncatts(self, attdict)
    
    def subsetVariables(self, varkeys, inplace = False):
        if not 'TFLAG' in varkeys:
            varkeys = ['TFLAG'] + list(varkeys)
        varkeylen = 16
        varliststr = getattr(self, 'VAR-LIST')
        varlist = list(varliststr[0+i:varkeylen+i].strip() for i in range(0, len(varliststr), varkeylen))
        outf = PseudoNetCDFFile.subsetVariables(self, varkeys, inplace = inplace)
        sliceo = [vi for vi, varkey in enumerate(varlist) if varkey in varkeys]
        newvarlist = [varkey for varkey in varlist if varkey in varkeys]
        outf = outf.sliceDimensions(VAR = sliceo)
        outf.NVAR = len(outf.dimensions['VAR'])
        setattr(outf, 'VAR-LIST', ''.join([vk.ljust(16) for vk in newvarlist]))
        return outf
    
    def sliceDimensions(self, *args, **kwds):
        import numpy as np
        outf = PseudoNetCDFFile.sliceDimensions(self, *args, **kwds)
        if 'LAY' in kwds:
            outf.VGLVLS[kwds['LAY']]
        if 'COL' in kwds and 'COL' in outf.dimensions:
            outf.NCOLS = len(outf.dimensions['COL'])
            outf.XORIG += np.r_[kwds['COL']][0] * outf.XCELL
        if 'ROW' in kwds and 'COL' in outf.dimensions:
            outf.NROWS = len(outf.dimensions['ROW'])
            outf.YORIG += np.r_[kwds['ROW']][0] * outf.YCELL
        if 'TSTEP' in kwds:
            import datetime
            times = np.atleast_1d(self.getTimes()[kwds['TSTEP']])
            outf.SDATE = int(times[0].strftime('%Y%j'))
            outf.STIME = int(times[0].strftime('%H%M%S'))
            if times.size > 1:
                dt = np.diff(times)
                if not (dt[0] == dt).all():
                    warn('New time is unstructured')
                outf.TSTEP = int((datetime.datetime(1900, 1, 1, 0) + dt[0]).strftime('%H%M%S'))
        return outf

class ioapi(Dataset, ioapi_base):
    def __getattribute__(self, *args, **kwds):
        return Dataset.__getattribute__(self, *args, **kwds)

    @classmethod
    def isMine(cls, *args, **kwds):
        try:
            f = Dataset(*args, **kwds)
            for dk in ['TSTEP', 'VAR', 'DATE-TIME']:
                assert(dk in f.dimensions)
            attrlist = f.ncattrs()
            for pk in ['XORIG', 'XCELL', 'YCELL', 'YORIG', 'SDATE', 'STIME']:
                assert(pk in attrlist)
            return True
        except:
            return False
