import PseudoNetCDF as pnc
from PseudoNetCDF.pncwarn import warn
import os
import io
import re
from datetime import datetime
from collections import OrderedDict
import numpy as np
from ._ioapi import ioapi_base

_prjp = ('GDTYP', 'P_ALP', 'P_BET', 'P_GAM', 'XCENT', 'YCENT')
_grdp = (
    'PRJNAME', 'XORIG', 'YORIG', 'XCELL', 'YCELL', 'NCOLS', 'NROWS', 'NTHIK'
)


class griddesc(ioapi_base):
    """
    griddesc is designed to read IOAPI griddesc files and create
    a file with basic CF metadata. An example format of GRIDDESC
    is show below with two grids on one projection.

' '
'LCC'
  2        33.000        45.000       -97.000       -97.000        40.000
' '
'SE52BENCH'
'LCC'    792000.000  -1080000.000     12000.000     12000.000 100  80   1
'36US3'
'LCC'  -2952000.000  -2772000.000     36000.000     36000.000 172  148   1
' '
"""
    def __init__(
        self, path, GDNAM=None, VGLVLS=(1., 0.), VGTOP=5000., FTYPE=1, VGTYP=7,
        SDATE=-635, STIME=0, TSTEP=0, var_kwds=None, nsteps=1, withcf=True,
        **prop_kw
    ):
        """
        Arguments
        ---------
        path : str or file-like
            If path is has a read method, it will be used directly.
            If path is a path to a file on disk, it will be read.
            If none of these, treat the path as text content.
            If path is None, load a griddesc using properties provided.
        GDNAM : str
            Name the grid to be used. If not provided, the first will be used
        VGLVLS : tuple
            Iterable of layer edge (k+1) values for the vertical grid.
        VGTOP : float
            Top of vertical grid.
        VGTYP : int
            Determines the units of VGLVLS and VGTOP. (7: WRF sigma-p,
            6: meters asl, for more details see IOAPI documentation)
        SDATE : int
            Starting julian date (YYYYJJJ) or -635
        STIME : int
            Starting time as HHMMSS
        TSTEP : int
            IOAPI time step in HHMMSS
        FTYPE : int
            1 for gridded; 2 for boundary
        var_kwds: tuple, list, or dict
            See setdefvars; Defaults to {'DUMMY': {'units': 'unknown'}}
        nsteps : int
            Number of time steps to use. Should be coupled with SDATE, STIME,
            and TSTEP. nsteps > 1 is not compatible with time independent
            (TSTEP=0) or SDATE=-635
        withcf : bool
            If true, then CF compatible variables that describe dimensions
            and time are added.
        prop_kw: mappable
            Can provide additional IOAPI (or other) properties for the file.
        """
        now = datetime.now()
        prj = {}
        grd = OrderedDict()
        if hasattr(path, 'read'):
            gdf = path
            path = '<inline>'
        elif path is None:
            required_props = _prjp[:] + _grdp[1:] + ('GDNAM',)
            grid_kw = prop_kw.copy()
            if GDNAM is None:
                grid_kw['GDNAM'] = 'UNKNOWN'
            else:
                grid_kw['GDNAM'] = GDNAM
            grid_kw.setdefault('NTHIK', 1)
            grid_kw.setdefault('PRJNAME', 'UNKNOWN')
            for rpk in required_props:
                if rpk not in grid_kw:
                    if path is None and 'GRIDDESC' in os.environ:
                        path = os.environ['GRIDDESC']
                        gdf = open(path, 'r')
                        break
                    raise KeyError(
                        f'{rpk} not found.'
                        + f'If path is None, {required_props} are required'
                    )
            else:
                gdf = io.StringIO("""
' '
'{PRJNAME}'
{GDTYP} {P_ALP} {P_BET} {P_GAM} {XCENT} {YCENT}
' '
'{GDNAM}'
'{PRJNAME}' {XORIG} {YORIG} {XCELL} {YCELL} {NCOLS} {NROWS} {NTHIK}
' '""".format(**grid_kw))
            path = '<inline>'
        elif os.path.exists(path):
            gdf = open(path, 'r')
        else:
            gdf = io.StringIO(path)
            path = '<inline>'

        gdtxt = gdf.read().replace(',', ' ').strip()
        # Fortran allows exponential notation to use D or d
        # instead of E or e, while Python does not.
        dble = re.compile('([\d.])[Dd]([\d])')
        gdtxt = dble.sub(r'\1e\2', gdtxt)
        gdlines = gdtxt.split('\n')
        gdlines = [line.strip() for line in gdlines]
        # Remove comments that start with !
        gdlines = [line.split('!')[0].strip() for line in gdlines]
        # IOAPI does not verify first line, simply discards.
        # dscgrid.f#L153
        # assert (gdlines[0].replace(' ', '') == "''")
        assert (gdlines[-1].replace(' ', '') == "''")
        i = 0
        blanks = []
        while i < len(gdlines):
            line = gdlines[i]
            if line.strip() in ("' '", "''", ""):
                blanks.append(i)
            else:
                i += 1
                parts = [eval(p) for p in gdlines[i].split()]
                key = eval(line)
                if len(blanks) == 1:
                    prj[key] = dict(zip(_prjp, parts))
                elif len(blanks) == 2:
                    grd[key] = dict(zip(_grdp, parts))
                else:
                    pass
            i += 1
        self._prj = prj
        self._grd = grd
        self.FTYPE = FTYPE
        self.VGTYP = VGTYP
        self.VGTOP = np.float32(VGTOP)
        self.VGLVLS = np.array(VGLVLS, dtype='f')
        self.NLAYS = len(VGLVLS) - 1
        self.UPNAM = 'GRIDDESC'
        self.EXEC_ID = 'GRIDDESC'
        self.FILEDESC = 'GRIDDESC'
        self.HISTORY = 'made from ' + path
        self.CDATE = int(now.strftime('%Y%j'))
        self.CTIME = int(now.strftime('%H%M%S'))
        self.WDATE = self.CDATE
        self.WTIME = self.CTIME
        self.SDATE = SDATE
        self.STIME = STIME
        self.TSTEP = TSTEP
        self._synthvars = None
        # user supplied properties should overwrite
        # defaults above
        for pk, pv in prop_kw.items():
            self.setncattr(pk, pv)

        self.setdefvars(var_kwds)
        self.addtime(nsteps)
        self.setgrid(key=GDNAM, withcf=withcf)
        self.updatemeta()
        self.getVarlist()

    def setdefvars(self, var_kwds):
        """
        Arguments
        ---------
        var_kwds: list, tuple, or dictionary
            Variables with units or just a list of variables
        """

        if var_kwds is None:
            var_kwds = {'DUMMY': {'units': 'unknown'}}
        elif isinstance(var_kwds, (tuple, list)):
            var_kwds = {k: 'unknown' for k in var_kwds}
        elif hasattr(var_kwds, 'items'):
            pass
        else:
            raise TypeError('var_kwds should be a list, tuple, or dictionary')

        self._var_kwds = {}
        for k, v in var_kwds.items():
            if isinstance(v, str):
                self._var_kwds[k] = {'units': v.ljust(16)}
            elif isinstance(v, dict):
                self._var_kwds[k] = {_k: _v for _k, _v in v.items()}
            else:
                raise TypeError(
                    'var_kwds should have either str values (units) or a dict'
                    + ' of all properties'
                )

    def setgrid(self, key=None, withcf=True):
        """
        Remakes the file to use the grid specified by key

        Arguments
        ---------
        key : str
            GDNAM to set the grid of the file
        withcf : bool
            Passed to adddims

        Returns
        -------
        None
        """
        if key is None:
            if getattr(self, 'GDNAM', None) is None:
                self.GDNAM = list(self._grd)[0]

            key = self.GDNAM.strip()
        grd = self._grd[key]
        prj = self._prj[grd['PRJNAME'].strip()]
        self.setncattr('GDNAM', key)
        self.setncatts(prj)
        self.setncatts(grd)
        self.adddims(withcf=withcf)

    def addtime(self, nsteps=1):
        """
        Adds TFLAG variable

        Arguments
        ---------
        nsteps : int
            Number of time steps to add

        Returns
        -------
        None
        """
        self.createDimension('TSTEP', nsteps).setunlimited(True)
        self.createDimension('DATE-TIME', 2)
        self.NVARS = len(self._var_kwds)
        self.createDimension('VAR', self.NVARS)
        try:
            self.updatetflag(overwrite=True)
        except Exception as e:
            warn(str(e))
            tflag = self.createVariable(
                'TFLAG', 'i',
                ('TSTEP', 'VAR', 'DATE-TIME'),
                units='<YYYYJJJ,HHMMSS>',
                long_name='TFLAG'.ljust(16),
                var_desc='TFLAG'.ljust(80)
            )
            tflag[:, :, 0] = self.SDATE
            tflag[:, :, 1] = 0

    def adddims(self, withcf=True):
        """
        Add spatial dimensions

        Arguments
        ---------
        withcf : bool
            If true, add Climate and Forecasting convention variables

        Returns
        -------
        None
        """
        for k in 'LAY ROW COL PERIM'.split():
            if k in self.dimensions:
                del self.dimensions[k]
        self.createDimension('LAY', self.NLAYS)
        if self.FTYPE == 1:
            self.createDimension('ROW', self.NROWS)
            self.createDimension('COL', self.NCOLS)
            dims = ('TSTEP', 'LAY', 'ROW', 'COL')
        elif self.FTYPE == 2:
            perim = (self.NCOLS * 2 + self.NROWS * 2 + self.NTHIK * 4)
            self.createDimension('PERIM', perim * self.NTHIK)
            dims = ('TSTEP', 'LAY', 'PERIM')
        else:
            raise ValueError(
                'Only supports FTYPE 1 or 2; received ' +
                str(self.FTYPE)
            )

        for vark, varkw in self._var_kwds.items():
            self.createVariable(
                vark, 'f', dims,
                long_name=vark.ljust(16),
                var_desc=vark.ljust(80),
                **varkw
            )

        if self._synthvars is None:
            oldkeys = set(self.variables)
        else:
            for key in self._synthvars:
                if key in self.variables:
                    del self.variables[key]
        if withcf:
            pnc.conventions.ioapi.add_cf_from_ioapi(self)
        if self._synthvars is None:
            self._synthvars = set(self.variables).difference(oldkeys)


if __name__ == '__main__':
    f = griddesc('GRIDDESC', '36US3')
