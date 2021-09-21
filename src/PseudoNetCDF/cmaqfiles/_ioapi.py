from PseudoNetCDF.pncwarn import warn
from ..core._files import PseudoNetCDFFile, netcdf
from collections import OrderedDict
from PseudoNetCDF._getwriter import registerwriter

import numpy as np
import datetime
today = datetime.datetime.today()
# Assuming some fo the most common
# options
_i = _ioapi_defaults = OrderedDict()
_i['IOAPI_VERSION'] = "N/A".ljust(80)
_i['EXEC_ID'] = "????????????????".ljust(80)
_i['FTYPE'] = 1
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


def ioapi_sort_meta(infile):
    mydimensions = infile.dimensions.copy()
    outvars = getattr(infile, 'VAR-LIST', '').split()
    allvars = outvars + \
        [k for k in list(infile.variables)
         if k not in outvars and k != 'TFLAG']
    infile.NVARS = len(outvars)
    infile.dimensions = OrderedDict()
    for dk in 'TSTEP DATE-TIME LAY VAR ROW COL'.split():
        dv = mydimensions[dk]
        if dk == 'VAR':
            dl = infile.NVARS
        else:
            dl = len(dv)
        ndv = infile.createDimension(dk, dl)
        ndv.setunlimited(dv.isunlimited())
    myvariables = infile.variables
    infile.variables = OrderedDict()
    for vk in ['TFLAG'] + allvars:
        infile.variables[vk] = myvariables[vk]


class ioapi_base(PseudoNetCDFFile):
    @classmethod
    def isMine(self, path, *args, **kwds):
        return False

    def __init__(self, *args, **kwds):
        super().__init__(self, *args, **kwds)
        self.setCoords(['TFLAG'])

    def __setattr__(self, k, v):
        """
        IOAPI has expected data types, which largely conform to the NetCDF3
        classic data model. For example, integers are int32 and floats are
        float64. However, VGTOP and VGLVLS are exceptions where the data type
        is expected to be float32. Further, any character variables will have
        one of three lengths: 16, 80, 80 * 60.
        """
        if k == 'VGTOP':
            v = np.float32(v)
        elif k == 'VGLVLS':
            v = np.asarray(v, dtype='f')
        elif k in ('FILEDESC', 'HISTORY'):
            if len(v) > 4800:
                warn(f'Truncating {k} 4800')
                v = v[:4800]
            v = v.ljust(4800)
        elif k in ('GDNAM', 'UPNAM'):
            if len(v) > 16:
                warn(f'Truncating {k} to 16')
                v = v[:16]
            v = v.ljust(16)
        elif k in ('EXEC_ID',):
            if len(v) > 80:
                warn(f'Truncating {k} to 80')
                v = v[:80]
            v = v.ljust(80)

        super().__setattr__(k, v)

    def audit_meta(self, fail='warn'):
        """
        Audit the IOAPI metadata. Checks existence of properties, dimensions
        internal consistency of properties and dimensions, and character
        lengths

        Arguments
        ---------
        fail : bool
           If fail, 'warn' or 'error' or 'ignore'

        Returns
        -------
        passing, audit, var_audit:
            audit stats and audit checklist dictionaries

        """
        audit = {}
        dimlens = {}
        for dk, dv in self.dimensions.items():
            dimlens[dk] = len(dv)

        for dk in ['TSTEP', 'DATE-TIME', 'PERIM', 'LAY', 'VAR', 'ROW', 'COL']:
            audit[f'has_{dk}'] = dk in dimlens
            if dk in ('TSTEP', 'DATE-TIME', 'PERIM'):
                continue
            audit[f'has_N{dk}S'] = hasattr(self, f'N{dk}S')
            if audit[f'has_N{dk}S']:
                audit[dk] = getattr(self, f'N{dk}S', 0) == dimlens[dk]

        if audit['has_ROW'] and audit['has_COL'] and not audit['has_PERIM']:
            del audit['has_PERIM']
            hdims = ('ROW', 'COL')
        elif audit['has_PERIM'] and not (audit['has_COL'] or audit['has_ROW']):
            del audit['has_ROW']
            del audit['has_COL']
            hdims = ('PERIM',)

        varliststr = getattr(self, 'VAR-LIST')
        audit['VAR-LIST-LEN'] = (self.NVARS * 16) == len(varliststr)
        checkvarliststr = self.getVarlist(update=False, retval='str')
        checkvarlist = self.getVarlist(update=False, retval='list')
        asisvarlist = self.getVarlist(update=False, prune=False, retval='list')
        for varkey in asisvarlist:
            audit[f'has_{varkey}'] = varkey in self.variables

        audit['VAR-LIST'] = varliststr == checkvarliststr
        var_audits = {}
        proplens = {
            'units': 16,
            'long_name': 16,
            'var_desc': 80,
        }
        for vk, var in self.variables.items():
            var_audit = var_audits[vk] = {}
            vardims = var.dimensions
            varshape = var.shape
            if vk == 'TFLAG':
                var_audit['right_dims'] = (
                    vardims == ('TSTEP', 'VAR', 'DATE-TIME')
                )
            else:
                var_audit['right_dims'] = (
                    vardims == ('TSTEP', 'LAY') + hdims
                )
            var_audit['ndims'] = len(varshape) == len(vardims)
            for dk, vn in zip(vardims, varshape):
                var_audit[f'shape_{dk}'] = vn == dimlens[dk]

            for pk, pcl in proplens.items():
                var_audit[f'has_{pk}'] = hasattr(var, pk)
                var_audit[f'len_{pk}'] = len(getattr(var, pk, '')) == pcl

            var_audit['SUMMARY'] = all(var_audit.values())
            if vk in checkvarlist:
                var_audit['INCLUDED'] = True
                audit[f'var_{vk}'] = var_audit['SUMMARY']
            else:
                var_audit['INCLUDED'] = False

        for pk, pcv in _ioapi_defaults.items():
            hasp = hasattr(self, pk)
            audit[f'has_{pk}'] = hasp
            if hasp:
                typep = isinstance(getattr(self, pk), type(pcv))
            else:
                typep = False
            audit[f'type_{pk}'] = typep

        proplens = {
            'FILEDESC': 4800,
            'HISTORY': 4800
        }
        for pk, pcl in proplens.items():
            audit[f'has_{pk}'] = hasattr(self, pk)
            audit[f'len_{pk}'] = len(getattr(self, pk, '')) <= pcl

        hastflag = 'TFLAG' in self.variables
        if hastflag:
            t0 = self.variables['TFLAG'][0, 0]
            audit['SDATE_TFLAG'] = (t0[0] == self.SDATE) or self.SDATE == -635
            audit['STIME_TFLAG'] = t0[1] == self.STIME

        audit['has_TFLAG'] = hastflag
        passing = audit['SUMMARY'] = all(list(audit.values()))
        if not passing:
            if fail == 'ignore':
                pass
            elif fail == 'warn':
                warn('Failed meta data audit review settings.')
            else:
                failed = str({k: v for k, v in audit.items() if not v})
                var_failed = str({
                    vk: {k: v for k, v in vv.items() if not v}
                    for vk, vv in var_audits.items()
                })
                raise ValueError(f"""Failed audit.
File failures: {failed}
Varable failures: {var_failed}
""")
        return passing, audit, var_audits

    def _updatetime(self, write=True, create=False):
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

    def createVariable(self, name, type='f', dimensions=None,
                       fill_value=None, **properties):
        """
        Wrapper on PseudoNetCDF.createVariable that updates VAR-LIST,
        NVARS, VAR, and TFLAG. Also adds long_name, var_desc, and units
        if not already in properties. long_name and var_desc default to
        name, while units defaults to unknown. These properties will
        also be adjusted to expected lengths (16, 80, 16).

        See also
        --------
        see PseudoNetCDFFile.createVariable
        """
        from copy import copy

        ftype = getattr(self, 'FTYPE', 1)
        hasperim = 'PERIM' in self.dimensions
        if dimensions is None:
            dimensions = ('TSTEP', 'LAY', 'ROW', 'COL')
            if (ftype == 2) and hasperim:
                dimensions = ('TSTEP', 'LAY', 'PERIM')

        properties = copy(properties)
        properties.setdefault('long_name', name.ljust(16))
        properties.setdefault('var_desc', name.ljust(80))
        properties.setdefault('units', 'unknown'.ljust(16))

        for pk, pv in list(properties.items()):
            if pk in ('long_name', 'units'):
                properties[pk] = pv.ljust(16)
            elif pk in ('var_desc',):
                properties[pk] = pv.ljust(80)

        if name == 'TFLAG':
            fill_value = None

        out = PseudoNetCDFFile.createVariable(
            self, name=name, type=type, dimensions=dimensions,
            fill_value=fill_value, **properties)

        self._add2Varlist([name])

        return out

    def copy(self, props=True, dimensions=True, variables=True, data=True):
        out = PseudoNetCDFFile.copy(
            self, props=props, dimensions=dimensions, variables=False,
            data=False
        )
        if variables:
            for vk, vv in self.variables.items():
                if not vk.endswith('TFLAG'):
                    PseudoNetCDFFile.copyVariable(
                        out, vv, key=vk, withdata=data
                    )

        if props and dimensions:
            out.updatetflag()

        return out

    def copyVariable(self, var, key=None, dtype=None, dimensions=None,
                     fill_value=None, withdata=True):
        """
        Wrapper on PseudoNetCDF.copyVariable that updates VAR-LIST,
        NVARS, VAR, and TFLAG

        See also
        --------
        see PseudoNetCDFFile.copyVariable
        """
        if key is None:
            for propk in ['_name', 'name', 'standard_name', 'long_name']:
                if hasattr(var, propk):
                    key = getattr(var, propk)
                    break
            else:
                raise AttributeError(
                    'varkey must be supplied because var has no name, ' +
                    'standard_name or long_name')

        outvar = PseudoNetCDFFile.copyVariable(
            self, var, key=key, dtype=dtype, dimensions=dimensions,
            fill_value=fill_value, withdata=withdata)

        # The long_name should be the same as the copied key
        outvar.long_name = key.ljust(16)
        # The var_desc and units should default to the copied variable
        if not hasattr(outvar, 'var_desc'):
            outvar.var_desc = key.ljust(80)
        if not hasattr(outvar, 'units'):
            outvar.units = 'unknown'.ljust(16)

        self._add2Varlist([key])
        return outvar

    def mask(self, *args, **kwds):
        """
        Wrapper on PseudoNetCDFFile.subsetVariables that updates VAR-LIST,
        NVARS, VAR, and TFLAG

        See also
        --------
        see PseudoNetCDFFile.mask
        """
        outf = PseudoNetCDFFile.mask(self, *args, **kwds)
        PseudoNetCDFFile.copyVariable(
            outf, self.variables['TFLAG'], key='TFLAG'
        )
        outf.updatemeta()
        return outf

    def subsetVariables(
        self, varkeys, inplace=False, exclude=False, keepcoords=True
    ):
        """
        Wrapper on PseudoNetCDFFile.subsetVariables that updates VAR-LIST,
        NVARS, VAR, and TFLAG

        See also
        --------
        see PseudoNetCDFFile.subsetVariables
        """
        varlist = self.getVarlist(update=False)
        newvarlist = [
            varkey for varkey in varlist
            if (
                (varkey in varkeys) != exclude and
                varkey in self.variables and
                not varkey.endswith('TFLAG')
            )
        ]
        outf = self.copy(props=True, dimensions=False, variables=False)
        for dk, dv in self.dimensions.items():
            if dk == 'VAR':
                outf.copyDimension(dv, key=dk, dimlen=len(newvarlist))
            else:
                outf.copyDimension(dv, key=dk)

        for vk in newvarlist:
            PseudoNetCDFFile.copyVariable(
                outf,
                self.variables[vk],
                key=vk,
                withdata=False
            )

        for vk in newvarlist:
            outf.variables[vk][:] = self.variables[vk][:]

        setattr(outf, 'VAR-LIST', '')
        outf._add2Varlist(newvarlist)
        outf.updatemeta()
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
        isarray = {
            dk: not np.isscalar(dv) and not isinstance(dv, slice)
            for dk, dv in dimslices.items()
        }
        # anyisarray = np.sum(list(isarray.values())) > 1

        # Check if COL or ROW was used
        hascol = 'COL' in dimslices
        hasrow = 'ROW' in dimslices
        deleterowcol = False
        if hascol and hasrow:
            if isarray['ROW'] and isarray['COL']:
                newdims = kwds.get('newdims', ('POINTS',))
                if 'ROW' not in newdims and 'COL' not in newdims:
                    deleterowcol = True

        # If lay was subset, subset VGLVLS too
        if 'LAY' in kwds:
            nlvls = outf.VGLVLS.size
            lidx = np.array(
                np.arange(outf.VGLVLS.size - 1)[kwds['LAY']], ndmin=1
            )
            tmpvglvls = outf.VGLVLS[lidx]
            if lidx[-1] < (nlvls - 1):
                try:
                    endl = outf.VGLVLS[lidx[-1] + 1]
                    tmpvglvls = np.append(tmpvglvls, endl)
                    outf.VGLVLS = tmpvglvls
                except Exception:
                    warn('VGLVLS could not be diagnosed; update manually')
        # If subsetting replaces ('ROW', 'COL') ... for example with ('PERIM',)
        # remove the dimensions
        if deleterowcol:
            del outf.dimensions['COL']
            del outf.dimensions['ROW']
        else:
            # Update origins
            if 'COL' in kwds and 'COL' in outf.dimensions:
                ncol = len(self.dimensions['COL'])
                outf.XORIG += np.arange(ncol)[kwds['COL']].take(0) * outf.XCELL
            if 'ROW' in kwds and 'ROW' in outf.dimensions:
                nrow = len(self.dimensions['ROW'])
                outf.YORIG += np.arange(nrow)[kwds['ROW']].take(0) * outf.YCELL

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
                outf.TSTEP = int(
                    (datetime.datetime(1900, 1, 1, 0) +
                     dt[0]).strftime('%H%M%S'))

        outf.updatemeta()
        return outf

    def interpSigma(self, vglvls, vgtop=None, interptype='linear',
                    extrapolate=False, fill_value='extrapolate',
                    verbose=0):
        """
        Parameters
        ----------
        vglvls : iterable
            the new vglvls (edges)
        vgtop : scalar
            Converting to new vgtop
        interptype : string
             'linear' uses a linear interpolation
             'conserve' uses a mass conserving interpolation
        extrapolate : boolean
            allow extrapolation beyond bounds with linear, default False
        fill_value : boolean
            set fill value (e.g, nan) to prevent extrapolation or edge
            continuation

        Returns
        -------
        outf : ioapi_base
            PseudoNetCDFFile with all variables interpolated

        Notes
        -----
        When extrapolate is false, the edge values are used for points beyond
        the inputs.
        """
        vglvls = np.asarray(vglvls)
        # grab input sigma coordinates
        myvglvls = self.VGLVLS

        # If needed, recalculate this files SIGMA
        if vgtop is not None and not vgtop == self.VGTOP:
            dp0 = 101325. - self.VGTOP
            dp1 = 101325. - vgtop
            myvglvls = (myvglvls * dp0 + self.VGTOP - vgtop) / dp1

        # Use midpoint for sigmas in inputs and outputs
        zs = (myvglvls[:-1] + myvglvls[1:]) / 2.
        nzs = (vglvls[:-1] + vglvls[1:]) / 2.

        if interptype == 'linear':
            from ..coordutil import getinterpweights
            weights = getinterpweights(
                zs, nzs, kind=interptype, fill_value=fill_value,
                extrapolate=extrapolate)
            # Create a function for interpolation

            def interpsigma(data):
                newdata = (weights * data[:, None]).sum(0)
                return newdata

        elif interptype == 'conserve':
            from ..coordutil import sigma2coeff
            # Calculate a weighting matrix using mass conserving
            # methods
            coeff = sigma2coeff(myvglvls, vglvls)  # (Nold, Nnew)
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
            raise ValueError(
                'interptype only implemented for "linear" and "conserve"')

        # Apply function on LAY
        outf = self.applyAlongDimensions(LAY=interpsigma, verbose=verbose)

        # Ensure vglvls is a simple array
        outf.VGLVLS = vglvls.view(np.ndarray).astype('f')
        outf.NLAYS = len(outf.VGLVLS) - 1
        outf.updatemeta()
        return outf

    def _add2Varlist(self, varkeys):
        varliststr = getattr(self, 'VAR-LIST', '')
        keys = [k for k in varliststr.split() if k in self.variables]
        newkeys = set(varkeys).difference(keys + ['ETFLAG', 'TFLAG'])
        for varkey in varkeys:
            if varkey in newkeys:
                varliststr += varkey.ljust(16)
                keys.append(varkey)
        setattr(self, 'NVARS', len(keys))
        setattr(self, 'VAR-LIST', varliststr)
        self._updatetime()
        return keys

    def getVarlist(self, update=True, prune=True, retval='list'):
        """
        Returns
        -------
        update : boolean
            update files attributes to be consistent
        prune : bool
            If True, remove variables that are missing or have
            names longer than 16 characters
        retval : str
            Return formatted string 'str' or list otherwise

        Notes
        -----
        If VAR-LIST does not exist, it is added assuming all variables
        with dimensions ('TSTEP', 'LAY', ...) are variables
        """
        if getattr(self, 'VAR-LIST', ''):
            varliststr_old = ''
            varlist = [k for k, v in self.variables.items()
                       if v.dimensions[:2] == ('TSTEP', 'LAY')]
        else:
            varliststr_old = getattr(self, 'VAR-LIST')
            if len(varliststr_old) % 16 == 0:
                ivars = int(len(varliststr_old) // 16)
                varlist = [
                    varliststr_old[ivar * 16:ivar * 16 + 16].strip()
                    for ivar in range(ivars)
                ]
            else:
                varlist = varliststr_old.split()

        # To be in the VAR-LIST, the variable must be shorter than 16
        # and be present and have the typical IOAPI dimensions
        vardimchecks = {}
        for vk in varlist:
            if vk in self.variables:
                dims = tuple(self.variables[vk].dimensions)
            else:
                dims = ()
            check = (
                dims == ('TSTEP', 'LAY', 'ROW', 'COL')
                or dims == ('TSTEP', 'LAY', 'PERIM')
            )
            if check and len(vk) > 16:
                warn(
                    f'{vk} name is too long ({len(vk)});'
                    + ' not included in VAR-LIST'
                )
                check = False

            vardimchecks[vk] = check

        varlist = [
            vk for vk in varlist
            if prune and (
                vardimchecks[vk]
                and len(vk) <= 16
            )
        ]

        varliststr_new = ''.join([vk.ljust(16)[:16] for vk in varlist])
        if update and varliststr_new != varliststr_old:
            setattr(self, 'VAR-LIST', varliststr_new)
        if update and len(varlist) != self.NVARS:
            self.NVARS = len(varlist)

        newdimlen = max(self.NVARS, 1)
        if 'VAR' in self.dimensions:
            if newdimlen != len(self.dimensions['VAR']):
                try:
                    self.createDimension('VAR', newdimlen)
                except Exception:
                    pass
                # add updatetflag
        else:
            self.createDimension('VAR', newdimlen)

        if retval == 'str':
            return varliststr_new
        else:
            return varlist

    def updatetflag(self, overwrite=None, startdate=None, tstep=None):
        if overwrite is None:
            overwrite = (
                'TFLAG' not in self.variables or
                self.variables['TFLAG'].shape[1] != self.NVARS
            )

        if overwrite:
            if 'TFLAG' in self.variables:
                del self.variables['TFLAG']
            if startdate is not None:
                self.SDATE = int(startdate.strftime('%Y%j'))
                self.STIME = int(startdate.strftime('%H%M%S'))
            if tstep is not None:
                self.TSTEP = tstep

            times = self.getTimes()
            tvar = self.createVariable(
                'TFLAG', 'i', ('TSTEP', 'VAR', 'DATE-TIME'))
            tvar.units = '<YYYYDDD,HHMMSS>'.ljust(16)
            tvar.long_name = 'TFLAG'.ljust(16)
            tvar.var_desc = ("Timestep-valid flags:  (1) YYYYDDD or (2) " +
                             "HHMMSS                                ")

            yyyyjjj = np.array([int(t.strftime('%Y%j')) for t in times])
            hhmmss = np.array([int(t.strftime('%H%M%S')) for t in times])

            tvar[:, :, 0] = yyyyjjj[:, None].repeat(tvar.shape[1], 1)
            tvar[:, :, 1] = hhmmss[:, None].repeat(tvar.shape[1], 1)
            self.SDATE = tvar[0, 0, 0]
            self.STIME = tvar[0, 0, 1]
        else:
            if len(self.dimensions['VAR']) == 0:
                return

            try:
                times = self.getTimes()
            except Exception as e:
                warn('Times were incalculable: using epoch start\n' + str(e))
                times = np.array([datetime.datetime(1970, 1, 1)])

            if not hasattr(self, 'SDATE'):
                self.SDATE = int(times[0].strftime('%Y%j'))
            if not hasattr(self, 'STIME'):
                self.STIME = int(times[0].strftime('%H%M%S'))
            if not hasattr(self, 'TSTEP'):
                if times.size > 1:
                    dt = np.diff(times)
                    if not (dt[0] == dt).all():
                        warn('New time is unstructured')

                    tstep = int(
                        (
                            datetime.datetime(1900, 1, 1, 0) + dt.mean()
                        ).strftime('%H%M%S')
                    )
                    self.TSTEP = tstep
                else:
                    self.TSTEP = 10000

    def updatemeta(self, attdict={}, sortmeta=False):
        """
        Parameters
        ----------
        attdict : dictionary
            key value pairs to update meta data
        sortmeta : boolean
            sort meta data after update
        Returns
        -------
        None

        Notes
        -----
        Meta data not provided or present will be inferred or made up.
        (See _ioapi_defaults)
        """
        if attdict is None:
            attdict = {}
        attdict.update(_ioapi_defaults)
        for pk, pv in attdict.items():
            if not hasattr(self, pk):
                setattr(self, pk, pv)

        if 'TSTEP' in self.dimensions:
            td = self.dimensions['TSTEP']
            if not td.isunlimited():
                td.setunlimited(True)

        if 'DATE-TIME' not in self.dimensions:
            self.createDimension('DATE-TIME', 2)

        self.getVarlist(update=True)

        if 'LAY' in self.dimensions:
            self.NLAYS = len(self.dimensions['LAY'])
        if 'COL' in self.dimensions:
            self.NCOLS = len(self.dimensions['COL'])
        if 'ROW' in self.dimensions:
            self.NROWS = len(self.dimensions['ROW'])

        self._updatetime()
        self.updatetflag()
        if sortmeta:
            ioapi_sort_meta(self)

    def applyAlongDimensions(self, *args, **kwds):
        """
        Wrapper PseudoNetCDFFile.applyAlongDimensions that corrects ROW, COL,
        LAY and TIME meta-data according to the ioapi format

        Parameters
        ----------
        see PseudoNetCDFFile.applyAlongDimensions
        """
        outf = PseudoNetCDFFile.applyAlongDimensions(self, *args, **kwds)
        if 'LAY' in kwds:
            nlays = len(self.dimensions['LAY'])
            layf = PseudoNetCDFFile()
            layf.createDimension('lay', nlays)
            layf.createDimension('nv', 2)
            laym = layf.createVariable('lay', 'f', ('lay',))
            layb = layf.createVariable('lay_bounds', 'f', ('lay', 'nv'))
            layb[:, 0] = self.VGLVLS[:-1]
            layb[:, 1] = self.VGLVLS[1:]
            laym[:] = layb.mean(1)
            newlayf = layf.applyAlongDimensions(lay=kwds['LAY'])
            nlayb = newlayf.variables['lay_bounds']
            outf.VGLVLS = np.append(nlayb[:, 0], nlayb[:, 1]).view(np.ndarray)
        outf.updatemeta()
        return outf

    def eval(self, *args, **kwds):
        """
        Wrapper PseudoNetCDFFile.eval that corrects VAR-LIST
        and TFLAG meta-data according to the ioapi format

        Parameters
        ----------
        see PseudoNetCDFFile.eval
        """
        # oldkeys = set(self.variables)
        out = PseudoNetCDFFile.eval(self, *args, **kwds)
        outkeys = set(out.variables)
        # newkeys = outkeys.difference(oldkeys)
        # byekeys = oldkeys.difference(outkeys)
        out._add2Varlist(outkeys)
        out.updatemeta()
        return out

    def getMap(self, maptype='basemap_auto', **kwds):
        """
        Wrapper PseudoNetCDFFile.getMap that uses NCOLS, XCELL
        NROWS, and YCELL to calculate map bounds if basemap_auto

        Parameters
        ----------
        see PseudoNetCDFFile.getMap
        """
        if maptype.endswith('_auto'):
            if self.GDTYP == 1:
                lllon, lllat = self.XORIG, self.YORIG
                urlon = self.XORIG + self.NCOLS * self.XCELL
                urlat = self.YORIG + self.NROWS * self.YCELL
            else:
                lllon, lllat = self.xy2ll(0, 0)
                urlon, urlat = self.xy2ll(
                    self.NCOLS * self.XCELL, self.NROWS * self.YCELL)
            kwds.setdefault('llcrnrlon', lllon)
            kwds.setdefault('llcrnrlat', lllat)
            kwds.setdefault('urcrnrlon', urlon)
            kwds.setdefault('urcrnrlat', urlat)
            maptype = maptype[:-5]

        return PseudoNetCDFFile.getMap(self, maptype=maptype, **kwds)

    def plot(self, varkey, plottype=None, ax_kw=None,
             plot_kw=None, cbar_kw=None, map_kw=None, dimreduction='mean'):
        """
        Parameters
        ----------
        varkey : string
            the variable to plot
        plottype : string
            longitude-latitude, latitude-pressure, longitude-pressure,
            vertical-profile, time-longitude, time-latitude,
            time-pressure, default, last two dimensions in reverse order
        ax_kw : dictionary
            keywords for the axes to be created
        plot_kw : dictionary
            keywords for the plot (plot, scatter, or pcolormesh) to be
            created
        cbar_kw : dictionary or bool or None
            keywords for the colorbar; if True or None, use defaults.
            If False, do not create a colorbar
        map_kw : dictionary or bool or None
            keywords for the getMap routine, which is only used with
            map capable dimensions (ie, plottype='longitude-latitude')
            If True or None, use default configuration dict(countries=True,
            coastlines=True, states=False, counties=False). If False,
            do not draw a map.
        dimreduction : string or function
            dimensions not being used in the plot are removed
            using applyAlongDimensions(dimkey=dimreduction) where
            each dimenions.
        """

        import matplotlib.pyplot as plt
        from ..coordutil import getbounds

        if ax_kw is None:
            ax_kw = {}

        if plot_kw is None:
            plot_kw = {}

        if cbar_kw is None or cbar_kw is True:
            cbar_kw = {}

        if map_kw is None or map_kw is True:
            map_kw = {}

        apply2dim = {}
        var = self.variables[varkey]
        if plottype is None:
            vardims = var.dimensions
            if len(vardims) > 1:
                plottype = '-'.join(vardims[-2:][::-1])
            else:
                plottype = vardims[0] + '-profile'
        varunit = varkey
        if hasattr(var, 'units'):
            varunit += ' ' + var.units.strip()
        dimlens = dict([(dk, len(self.dimensions[dk]))
                        for dk in var.dimensions])
        dimpos = dict([(dk, di) for di, dk in enumerate(var.dimensions)])
        raw_xkey, raw_ykey = plottype.split('-')
        d2d = {'time': 'TSTEP', 'latitude': 'ROW',
               'longitude': 'COL', 'pressure': 'LAY'}
        xkey = d2d.get(raw_xkey, raw_xkey)
        ykey = d2d.get(raw_ykey, raw_ykey)
        if not ykey == 'profile':
            for dimkey in list(dimlens):
                if dimkey not in (xkey, ykey) and dimlens.get(dimkey, 1) > 1:
                    apply2dim[dimkey] = dimreduction

        if len(apply2dim) > 0:
            myf = self.applyAlongDimensions(**apply2dim)
            var = myf.variables[varkey]
            dimlens = dict([(dk, len(self.dimensions[dk]))
                            for dk in var.dimensions])
        else:
            myf = self
        if ykey in ('profile',):
            vaxi = var.dimensions.index(xkey)
            vsize = var.shape[vaxi]
            vals = np.rollaxis(var[:], vaxi).reshape(vsize, -1)
        else:
            vals = var[:].squeeze()

        if xkey == 'TSTEP':
            xm = myf.getTimes()
            dx = np.diff(xm)[-1]
            x = np.append(xm, xm[-1] + dx)
            x = plt.matplotlib.dates.date2num(x)
        else:
            x = getbounds(myf, xkey)

        ax = plt.gca(**ax_kw)
        if ykey in ('profile',):
            y = getbounds(myf, xkey)
            x1 = vals[:].min(1)
            xm = vals[:].mean(1)
            x2 = vals[:].max(1)
            ax.fill_betweenx(y=y, x1=x1, x2=x2, label=varkey + '(min, max)')
            ax.plot(xm, y, label=varkey, **plot_kw)
            ax.set_ylabel(xkey)
            ax.set_xlabel(varunit)
            return ax

        if ykey == 'TSTEP':
            ym = myf.getTimes()
            dy = np.diff(ym)[-1]
            y = np.append(ym, ym[-1] + dy)
            y = plt.matplotlib.dates.date2num(y)
        else:
            y = getbounds(myf, ykey)

        if dimpos[xkey] < dimpos[ykey]:
            vals = vals.T
        if xkey == 'TSTEP':
            ax.xaxis.set_major_formatter(
                plt.matplotlib.dates.AutoDateFormatter(
                    plt.matplotlib.dates.AutoDateLocator()))
        if ykey == 'TSTEP':
            ax.yaxis.set_major_formatter(
                plt.matplotlib.dates.AutoDateFormatter(
                    plt.matplotlib.dates.AutoDateLocator()))
        mappabledims = (
            plottype == 'longitude-latitude' or
            plottype == 'COL-ROW'
        )
        if mappabledims:
            if map_kw is not False:
                try:
                    coastlines = map_kw.pop('coastlines', True)
                    countries = map_kw.pop('countries', True)
                    states = map_kw.pop('states', False)
                    counties = map_kw.pop('counties', False)
                    bmap = myf.getMap(**map_kw)
                    if coastlines:
                        bmap.drawcoastlines(ax=ax)
                    if countries:
                        bmap.drawcountries(ax=ax)
                    if states:
                        bmap.drawstates(ax=ax)
                    if counties:
                        bmap.drawcounties(ax=ax)
                    x = np.arange(self.NCOLS + 1) * self.XCELL
                    y = np.arange(self.NROWS + 1) * self.YCELL
                    if self.GDTYP == 1:
                        x += self.XORIG
                        y += self.YORIG
                except Exception:
                    pass
        else:
            ax.set_xlabel(xkey)
            ax.set_ylabel(ykey)

        p = ax.pcolormesh(x, y, vals, **plot_kw)
        if cbar_kw is not False:
            cbar_kw.setdefault('label', varunit)
            ax.figure.colorbar(p, **cbar_kw)
        return ax

    slice = sliceDimensions
    apply = applyAlongDimensions


class ioapi(ioapi_base, netcdf):
    def __init__(self, *args, **kwds):
        netcdf.__init__(self, *args, **kwds)
        self.setCoords(['TFLAG'])

    def _newlike(self):
        if self.get_dest() is not None:
            outf = ioapi(**self.get_dest())
        elif isinstance(self, PseudoNetCDFFile):
            outt = ioapi_base
            outf = outt.__new__(outt)
        else:
            outf = PseudoNetCDFFile()
        outf.set_varopt(**self.get_varopt())
        outf._updatetime(write=True, create=True)
        return outf

    @classmethod
    def from_ncf(cls, infile):
        outf = ioapi_base()
        for pk in infile.ncattrs():
            pv = getattr(infile, pk)
            setattr(outf, pk, pv)

        for dk, dv in infile.dimensions.items():
            outf.copyDimension(dv, key=dk)

        for vk, vv in infile.variables.items():
            outf.copyVariable(vv, key=vk)

        return outf

    @property
    def _mode(self):
        return self.__dict__['_mode']

    def createVariable(self, *args, **kwds):
        return netcdf.createVariable(self, *args, **kwds)

    def createDimension(self, *args, **kwds):
        return netcdf.createDimension(self, *args, **kwds)

    @classmethod
    def isMine(cls, *args, **kwds):
        try:
            f = netcdf(*args, **kwds)
            for dk in ['TSTEP', 'VAR', 'DATE-TIME']:
                assert(dk in f.dimensions)
            attrlist = f.ncattrs()
            for pk in ['XORIG', 'XCELL', 'YCELL', 'YORIG', 'SDATE', 'STIME']:
                assert(pk in attrlist)
            return True
        except Exception:
            return False


def _fromdefaults(nfile, opts):
    for dk in opts:
        if dk in nfile.dimensions:
            return dk
    else:
        return opts[0]


def _tryset(obj, pk, pv, prefix='obj'):
    try:
        setattr(obj, pk, pv)
    except Exception as e:
        warn(
            '{}.{} not set {}; value {} was lost'.format(prefix, pk, e, pv)
        )


def ncf2ioapi(
    nfile, outpath, verbose=0, mode='w',
    tstepkeys=('TSTEP', 'Time', 'time', 't'),
    laykeys=('LAY', 'layer', 'level', 'layer47', 'layer72', 'z'),
    rowkeys=('ROW', 'latitude', 'south_north', 'lat', 'y'),
    colkeys=('COL', 'longitude', 'west_east', 'lon', 'x'),
    perimkeys=('PERIM', 'POINTS'),
    format='NETCDF3_CLASSIC',
    **props
):
    """
    Parameters
    ----------
    nfile : netcdf-like file
        Input file to translate to IOAPI meta data NETCDF file
    outpath : str
        Path for output file
    verbose : int
        Higher verbosity is more info
    mode : str
        Mode for opening file (w, ws, r+, a)
    tstepkeys : tuple of strings
        Strings to translate to TSTEP
    laykeys : tuple of strings
        Strings to translate to LAY
    rowkeys : tuple of strings
        Strings to translate to ROW
    colkeys : tuple of strings
        Strings to translate to COL
    perimkeys : tuple of strings
        Strings to translate to PERIM; only relevant when GDTYP=2
    format : string
        Format for output file either NETCDF3_CLASSIC, NETCDF3_64BIT_OFFSET
        NETCDF4_CLASSIC
    props : dict
        Mapping of properties to set

    Returns
    -------
    out : netCDF4.Dataset
    """
    time = _fromdefaults(nfile, tstepkeys)
    lay = _fromdefaults(nfile, laykeys)
    row = _fromdefaults(nfile, rowkeys)
    col = _fromdefaults(nfile, colkeys)
    perim = _fromdefaults(nfile, perimkeys)
    dim2dim = {}
    if time != 'TSTEP':
        dim2dim[time] = 'TSTEP'
    if lay != 'LAY':
        dim2dim[lay] = 'LAY'
    if row != 'ROW':
        dim2dim[row] = 'ROW'
    if col != 'COL':
        dim2dim[col] = 'COL'
    if col != 'COL':
        dim2dim[col] = 'COL'

    if verbose > 0:
        print('Dimension translation', dim2dim)

    fileprops = _ioapi_defaults.copy()
    for k in nfile.ncattrs():
        fileprops[k] = getattr(nfile, k)

    fileprops.update(**props)

    if fileprops['FTYPE'] == 2:
        indims = (time, lay, perim)
    else:
        indims = (time, lay, row, col)

    outdims = tuple(dim2dim.get(k, k) for k in indims)
    varkeys = [k for k, v in nfile.variables.items() if v.dimensions == indims]
    varlist = ''.join([k.ljust(16) for k in varkeys])

    fileprops['VAR-LIST'] = varlist
    nvar = fileprops['NVARS'] = len(varkeys)
    nthk = fileprops['NTHIK']
    nx = fileprops['NCOLS'] = len(nfile.dimensions[col])
    ny = fileprops['NROWS'] = len(nfile.dimensions[row])
    nz = fileprops['NLAYS'] = len(nfile.dimensions[lay])

    ofile = netcdf(outpath, format=format, mode=mode)

    for pk, pv in fileprops.items():
        _tryset(ofile, pk, pv, prefix='file')

    ofile.createDimension('TSTEP', None)
    ofile.createDimension('DATE-TIME', 2)
    ofile.createDimension('LAY', nz)
    ofile.createDimension('VAR', nvar)

    if fileprops['FTYPE'] == 2:
        ofile.createDimension('PERIM', nthk * (4 * nthk + 2 * (nx + ny)))
    elif fileprops['FTYPE'] == 1:
        ofile.createDimension('ROW', ny)
        ofile.createDimension('COL', nx)
    else:
        raise ValueError('FTYPE is unknown; must be 1 or 2')

    times = nfile.getTimes()
    tv = ofile.createVariable('TFLAG', 'i', ('TSTEP', 'VAR', 'DATE-TIME'))
    tv.units = '<YYYYJJJ,HHMMSS>'
    tv.long_name = 'TFLAG'.ljust(16)
    tv.var_desc = 'TFLAG'.ljust(80)
    yyyyjjj = np.array([t.strftime('%Y%j') for t in times], dtype='i')
    hhmmss = np.array([t.strftime('%H%M%S') for t in times], dtype='i')
    for k in varkeys:
        v = nfile.variables[k]
        ov = ofile.createVariable(k, v.dtype.char, outdims)
        ov.long_name = k.ljust(16)
        ov.var_desc = k.ljust(16)
        ov.units = 'unknown'.ljust(16)
        for pk in v.ncattrs():
            pv = getattr(v, pk)
            _tryset(ov, pk, pv, prefix=k)

    tv[:yyyyjjj.size, :, 0] = yyyyjjj[:, None].repeat(nvar, 1)
    tv[:yyyyjjj.size, :, 1] = hhmmss[:, None].repeat(nvar, 1)
    dt = (times[-1] - times[0]).total_seconds() / (len(times) - 1)
    tmpd = datetime.datetime(1900, 1, 1) + datetime.timedelta(seconds=dt)
    ofile.TSTEP = int(tmpd.strftime('%H%M%S'))
    ofile.SDATE = int(times[0].strftime('%Y%j'))
    ofile.STIME = int(times[0].strftime('%H%M%S'))
    # if (
    #     any([vk.startswith('lat') for vk in nfile.variables]) and
    #     any([vk.startswith('lon') for vk in nfile.variables])
    # ):
    #     ofile.GDTYP = 1
    #     projstr = nfile.getproj(projformat='proj4', withgrid=True)
    if not hasattr(ofile, 'VGLVLS'):
        ofile.VGTYP = -1
        ofile.VGLVLS = np.arange(nz + 1, dtype='f')
    for k in varkeys:
        v = nfile.variables[k]
        ov = ofile.variables[k]
        ov[:] = v

    return ofile


registerwriter('ioapi', ncf2ioapi)
