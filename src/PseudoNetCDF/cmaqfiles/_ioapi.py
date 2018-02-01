from PseudoNetCDF.pncwarn import warn
from netCDF4 import Dataset
from ..core._files import PseudoNetCDFFile
from collections import OrderedDict
import numpy as np
import datetime
today = datetime.datetime.today()
# Assuming some fo the most common
# options
_i = _ioapi_defaults = OrderedDict()
_i['IOAPI_VERSION'] = "N/A".ljust(80)
_i['EXEC_ID'] = "????????????????".ljust(80)
_i['FTYPE'] = -1
_i['CDATE'] = int(today.strftime('%Y%j'))
_i['CTIME'] = int(today.strftime('%H%M%S'))
_i['WDATE'] = int(today.strftime('%Y%j'))
_i['WTIME'] = int(today.strftime('%H%M%S'))
_i['NTHIK'] = 1
_i['GDTYP'] = 1
_i['P_ALP'] = 33.
_i['P_BET'] = 45.
_i['P_GAM'] = -97.
_i['XCENT'] = -97.
_i['YCENT'] = 40.
_i['VGTYP'] = 2
_i['VGTOP'] = np.float32(5000)
_i['GDNAM'] = "UNKNOWN         "
_i['UPNAM'] = "MAKEIOAPI       "
_i['FILEDESC'] = "".ljust(80)
_i['HISTORY'] = ""

class ioapi_base(PseudoNetCDFFile):
    def __getattributte__(self, *args, **kwds):
        return getattr(self, *args, **kwds)
    
    @classmethod
    def isMine(self, path):
        return False
    
    def _updatetime(self, write = True, create = False):
        from datetime import datetime
        t = datetime.now()
        try:
            if create:
                self.CDATE = int(t.strftime('%Y%j'))
                self.CTIME = int(t.strftime('%H%M%S'))
            if write:
                self.WDATE = int(t.strftime('%Y%j'))
                self.WTIME = int(t.strftime('%H%M%S'))
        except Exception as e:
            warn('Time could not be updated; ' + str(e))
        
    def setncatts(self, attdict):
        """
        Wrapper on PseudoNetCDF.setncatts that updates WDATE, and WTIME
        
        See also
        --------
        see PseudoNetCDFFile.setncatts
        """
        PseudoNetCDFFile.setncatts(self, attdict)
        self._updatetime()

    def createVariable(self, name, type, dimensions, fill_value = None, **properties):
        """
        Wrapper on PseudoNetCDF.createVariable that updates VAR-LIST,
        NVARS, VAR, and TFLAG
        
        See also
        --------
        see PseudoNetCDFFile.createVariable
        """
        out = PseudoNetCDFFile.createVariable(self, name = name, type = type, dimensions = dimensions, fill_value = fill_value, **properties)
        self._add2Varlist([name])
        return out
    
    def copyVariable(self, var, key = None, dtype = None, dimensions = None, fill_value = None, withdata = True):
        """
        Wrapper on PseudoNetCDF.copyVariable that updates VAR-LIST,
        NVARS, VAR, and TFLAG
        
        See also
        --------
        see PseudoNetCDFFile.copyVariable
        """
        outvar = PseudoNetCDFFile.copyVariable(self, var, key = key, dtype = dtype, dimensions = dimensions, fill_value = fill_value, withdata = withdata)
        if key is None:
            for propk in ['name', 'standard_name', 'long_name']:
                if hasattr(var, propk):
                    key = getattr(var, propk)
            else:
                raise AttributeError('varkey must be supplied because var has no name, standard_name or long_name')
        self._add2Varlist([key])
        return outvar
    
    def subsetVariables(self, varkeys, inplace = False, exclude = False):
        """
        Wrapper on PseudoNetCDFFile.subsetVariables that updates VAR-LIST,
        NVARS, VAR, and TFLAG
        
        See also
        --------
        see PseudoNetCDFFile.sliceDimensions
        """
        varlist = self.getVarlist(update = False)
        newvarlist = [varkey for varkey in varlist if (varkey in varkeys) != exclude]
        outf = PseudoNetCDFFile.subsetVariables(self, varkeys, inplace = inplace, exclude = exclude)
        if not 'TFLAG' in outf.variables and 'TFLAG' in self.variables:
            outf.copyVariable(self.variables['TFLAG'], key = 'TFLAG')
        sliceo = np.array([vi for vi, varkey in enumerate(varlist) if varkey in newvarlist])
        outf = outf.sliceDimensions(VAR = sliceo)
        setattr(outf, 'VAR-LIST', '')
        outf._add2Varlist(newvarlist)
        return outf
    
    def sliceDimensions(self, *args, **kwds):
        """
        Wrapper PseudoNetCDFFile.sliceDimensions that corrects ROW, COL,
        LAY and TIME meta-data according to the ioapi format
        
        Parameters
        ----------
        see PseudoNetCDFFile.sliceDimensions
        """
        # First slice as normal
        outf = PseudoNetCDFFile.sliceDimensions(self, *args, **kwds)
        # Copy slice keywords excluding newdims
        dimslices = kwds.copy()
        dimslices.pop('newdims', None)
        
        # Identify array indices and the need for fancy indexing
        isarray = {dk: not np.isscalar(dv) for dk, dv in dimslices.items()}
        anyisarray = np.sum(list(isarray.values())) > 1
        
        # Check if COL or ROW was used
        hascol = 'COL' in dimslices
        hasrow = 'ROW' in dimslices
        deleterowcol = False
        if hascol and hasrow:
            if isarray['ROW'] and isarray['COL']:
               deleterowcol = True
        
        # If lay was subset, subset VGLVLS too
        if 'LAY' in kwds:
            outf.VGLVLS = outf.VGLVLS[kwds['LAY']]
        
        # If subsetting replaces ('ROW', 'COL') ... for example with ('PERIM',)
        # remove the dimensions
        if deleterowcol:
            del outf.dimensions['COL']
            del outf.dimensions['ROW']
        else:
            # Update origins
            if 'COL' in kwds and 'COL' in outf.dimensions:
                outf.XORIG += np.r_[kwds['COL']][0] * outf.XCELL
            if 'ROW' in kwds and 'COL' in outf.dimensions:
                outf.YORIG += np.r_[kwds['ROW']][0] * outf.YCELL
        
        # Update TFLAG, SDATE, STIME and TSTEP
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
        
        outf.updatemeta()
        return outf
    
    def interpSigma(self, vglvls, vgtop = None, interptype = 'linear', extrapolate = False, fill_value = 'extrapolate'):
        """
        Parameters
        ----------
        self : the file to interpolate from must have VGLVLS
        vglvls : the new vglvls (edges)
        vgtop : Converting to new vgtop
        interptype : 'linear' or 'conserve'
             linear : uses a linear interpolation
             conserve : uses a mass conserving interpolation
        extrapolate : allow extrapolation beyond bounds with linear, default False
        fill_value : set fill value (e.g, nan) to prevent extrapolation or edge 
                     continuation
        
        Returns
        -------
        outf - ioapi_base PseudoNetCDFFile with al variables interpolated
        
        Notes
        -----
        When extrapolate is false, the edge values are used for points beyond
        the inputs.
        """
        # grab input sigma coordinates
        myvglvls = self.VGLVLS
        
        # If needed, recalculate this files SIGMA
        if not vgtop is None and not vgtop == self.VGTOP:
            dp0 = 101325. - self.VGTOP
            dp1 = 101325. - vgtop
            myvglvls = (myvglvls * dp0 + self.VGTOP - vgtop) / dp1
        
        # Use midpoint for sigmas in inputs and outputs
        zs = (myvglvls[:-1]+myvglvls[1:])/2.
        nzs = (vglvls[:-1]+vglvls[1:])/2.
        
        
        if interptype == 'linear':
            from ..coordutil import getinterpweights
            weights = getinterpweights(zs, nzs, kind = interptype, fill_value = fill_value, extrapolate = extrapolate)
            # Create a function for interpolation
            def interpsigma(data):
                newdata = (weights * data[:, None]).sum(0)
                return newdata
            
        elif interptype == 'conserve':
            from ..coordutil import sigma2coeff
            # Calculate a weighting matrix using mass conserving
            # methods
            coeff = sigma2coeff(myvglvls, vglvls) # (Nold, Nnew)
            # Calculate input mass fractions
            dp_in = -np.diff(myvglvls.astype('d'))[:, None]
            
            # Create a weighted mass fraction function
            fdp = dp_in * coeff
            
            # Create a precalculated normalizer
            ndp = fdp.sum(0)
            
            # Create a function for interpolation
            def interpsigma(data):
                nvals = (data[:, None] * fdp).sum(0) / ndp
                return nvals
        else:
            raise ValueError('interptype only implemented for "linear" and "conserve"')

        # Apply function on LAY
        outf = self.applyAlongDimensions(LAY = interpsigma)
        
        # Ensure vglvls is a simple array
        outf.VGLVLS = vglvls.view(np.ndarray).astype('f')
        outf.NLAYS = len(outf.VGLVLS) - 1
        outf.updatemeta()
        return outf
    
    def _add2Varlist(self, varkeys):
        varliststr = getattr(self, 'VAR-LIST', '')
        keys = [k for k in varliststr.split() if k in self.variables]
        newkeys = set(varkeys).difference(keys + ['TFLAG'])
        for varkey in varkeys:
            if varkey in newkeys:
                varliststr += varkey.ljust(16)
                keys.append(varkey)
        setattr(self, 'NVARS', len(keys))
        setattr(self, 'VAR-LIST', varliststr)
        self._updatetime()
        return keys
    
    def getVarlist(self, update = True):
        """
        Returns
        -------
        varlist : VAR-LIST split and stripped
        update  : update files attributes to be consistent
        
        Notes
        -----
        If VAR-LIST does not exist, it is added assuming all variables
        with dimensions ('TSTEP', 'LAY', ...) are variables
        """
        if not hasattr(self, 'VAR-LIST'):
            varliststr_old = ''
            varlist = ''.join([k.ljust(16) for k, v in self.variables.items() if v.dimensions[:2] == ('TSTEP', 'LAY')])
        else:
            varliststr_old = getattr(self, 'VAR-LIST')
            varlist = [vk for vk in varliststr_old.split() if vk in self.variables]
        varliststr_new = ''.join([vk.ljust(16) for vk in varlist])
        if update and varliststr_new != varliststr_old:
            setattr(self, 'VAR-LIST', varliststr_new)
        if update and len(varlist) != self.NVARS:
            self.NVARS = len(varlist)
            if 'VAR' in self.dimensions:
                if self.NVARS != len(self.dimensions['VAR']):
                    try: self.createDimension('VAR', self.NVARS)
                    except: pass
            else:
                self.createDimension('VAR', self.NVARS)

        return varlist
        
    def updatemeta(self, **attdict):
        """
        Parameters
        ----------
        attdict : key value pairs to update meta data
        
        Returns
        -------
        None
        
        Notes
        -----
        Meta data not provided or present will be inferred or made up. (See _ioapi_defaults)
        """
        attdict.update(_ioapi_defaults)
        for pk, pv in attdict.items():
            if not hasattr(self, pk):
                setattr(self, pk, pv)
        
        if 'TSTEP' in self.dimensions:
            td = self.dimensions['TSTEP']
            if not td.isunlimited():
                td.setunlimited(True)
        
        if 'LAY' in self.dimensions: self.NLAYS = len(self.dimensions['LAY'])
        if 'COL' in self.dimensions: self.NCOLS = len(self.dimensions['COL'])
        if 'ROW' in self.dimensions: self.NCOLS = len(self.dimensions['ROW'])
        
        self.getVarlist()
        
        if not 'DATE-TIME' in self.dimensions:
            self.createDimension('DATE-TIME', 2)
        
        self._updatetime()
        try:
            times = self.getTimes()
        except:
            return
        
        if not hasattr(self, 'SDATE'):
            self.SDATE = int(times[0].strftime('%Y%j'))
        if not hasattr(self, 'STIME'):
            self.STIME = int(times[0].strftime('%H%M%S'))
        if times.size > 1:
            dt = np.diff(times)
            if not (dt[0] == dt).all():
                warn('New time is unstructured')
            self.TSTEP = int((datetime.datetime(1900, 1, 1, 0) + dt[0]).strftime('%H%M%S'))
        if not 'TFLAG' in self.variables:
            dotflag = True
            tvar = self.createVariable('TFLAG', 'i', ('TSTEP', 'VAR', 'DATE-TIME'))
            tvar.units = '<YYYYDDD, HHMMSS>'.ljust(16)
            tvar.long_name = 'TFLAG'.ljust(16)
            tvar.var_desc = "Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS                                "
        else:
            tvar = self.variables['TFLAG']
            dotflag = ~((self.SDATE == tvar[0,0,0]) and (self.STIME == tvar[0,0,1]))
            
        
        if dotflag:
            yyyyjjj = np.array([int(t.strftime('%Y%j')) for t in times])
            hhmmss = np.array([int(t.strftime('%H%M%S')) for t in times])
            
            tvar[:, :, 0] = yyyyjjj[:, None]
            tvar[:, :, 1] = hhmmss[:, None]
        
    def applyAlongDimensions(self, *args, **kwds):
        outf = PseudoNetCDFFile.applyAlongDimensions(self, *args, **kwds)
        outf.updatemeta()
        return outf

class ioapi(Dataset, ioapi_base):
    def __getattribute__(self, *args, **kwds):
        return Dataset.__getattribute__(self, *args, **kwds)

    def _newlike(self):
        if isinstance(self, PseudoNetCDFFile):
            outt = ioapi_base
            outf = outt.__new__(outt)
        else:
            outf = PseudoNetCDFFile()
        outf._updatetime(write = True, create = True)
        return outf 
    
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
