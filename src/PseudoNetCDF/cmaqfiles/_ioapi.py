from netCDF4 import Dataset
from ..core._files import PseudoNetCDFFile

def _sigma2coeff(fromvglvls, tovglvls):
    """
    Calculate fraction of pressure from each layer in fromfile
    that is in each layer in tofile and return matrix
    """
    import numpy as np
    edges = np.interp(tovglvls[::-1], fromvglvls[::-1], np.arange(fromvglvls.size)[::-1])[::-1].repeat(2,0)[1:-1].reshape(-1, 2).astype('d')
    coeff = np.zeros((fromvglvls.size - 1, tovglvls.size - 1), dtype = 'd')
    for li, (b, t) in enumerate(edges):
        ll = np.floor(b).astype('i')
        ul = np.ceil(t).astype('i')
        for lay in range(ll, ul):
            bf = max(b - lay, 0)
            tf = min(t - lay, 1)
            myf = min(bf, tf)
            myf = tf - bf
            coeff[lay, li] = myf
            #if li == 12:
            #    print(li, b, t, ll, ul, lay, bf, tf, min(bf, tf))
    return coeff

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
        dimslices = kwds.copy()
        dimslices.pop('newdims', None)
        
        isarray = {dk: not np.isscalar(dv) for dk, dv in dimslices.items()}
        anyisarray = np.sum(list(isarray.values())) > 1
        hascol = 'COL' in dimslices
        hasrow = 'ROW' in dimslices
        deleterowcol = False
        if hascol and hasrow:
            if isarray['ROW'] and isarray['COL']:
               deleterowcol = True
        if 'LAY' in kwds:
            outf.VGLVLS[kwds['LAY']]
        if deleterowcol:
            del outf.dimensions['COL']
            del outf.dimensions['ROW']
        else:
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
    
    def interpSigma(self, vglvls, interptype = 'linear'):
        """
        self - the file to interpolate from must have VGLVLS
        vglvls - the new vglvls (edges)
        interptype - 'linear' or 'conserve'
             linear - uses a linear interpolation
             conserve - uses a mass conserving interpolation
        """
        import numpy as np
        from scipy.interpolate import interp1d
        # input sigma coordinates
        zs = (self.VGLVLS[:-1]+self.VGLVLS[1:])/2
        nzs = (vglvls[:-1]+vglvls[1:])/2
        if interptype == 'linear':
            # output sigma coordinates
            # identity matrix
            ident = np.identity(zs.size)
            # weight function
            weight_func = interp1d(zs, ident, 'linear')
            # weights
            weights = weight_func(nzs)
            def interpsigma(data):
                newdata = (weights * data[:, None]).sum(0)
                return newdata
        elif interptype == 'conserve':
            coeff = _sigma2coeff(self.VGLVLS, vglvls) # (Nold, Nnew)
            dp_in = -np.diff(self.VGLVLS.astype('d'))[:, None]
            #dp_out = -np.diff(VGLVLS.astype('d'))[:, None]
            fdp = dp_in * coeff
            ndp = fdp.sum(0)
            def interpsigma(data):
                nvals = (data[:, None] * fdp).sum(0) / ndp
                return nvals
        outf = self.applyAlongDimensions(LAY = interpsigma)
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
