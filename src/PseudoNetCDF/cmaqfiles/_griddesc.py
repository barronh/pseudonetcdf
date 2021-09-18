import PseudoNetCDF as pnc
import os
import io
from datetime import datetime
from collections import OrderedDict
import numpy as np
from ._ioapi import ioapi_base

_prjp = ['GDTYP', 'P_ALP', 'P_BET', 'P_GAM', 'XCENT', 'YCENT']
_grdp = [
    'PRJNAME', 'XORIG', 'YORIG', 'XCELL', 'YCELL', 'NCOLS', 'NROWS', 'NTHIK'
]


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
        withcf=True, **prop_kw
    ):
        """
        Arguments
        ---------
        path : str or file-like
            If path is has a read method, it will be used directly.
            If path is a path to a file on disk, it will be read.
            If none of these, it will try to treat the path as content.
        GDNAM : str
            Name the grid to be used. If not provided, the first will be used
        VGLVLS : tuple
            Iterable of layer edge (k+1) values for the vertical grid.
        VGTOP : float
            Top of vertical grid.
        VGTYP : int
            Determines the units of VGLVLS and VGTOP. (7: WRF sigma-p,
            6: meters asl, for more details see IOAPI documentation)
        FTYPE : int
            1 for gridded; 2 for boundary
        withcf : bool
            If true, then CF compatible variables that describe dimensions
            and time are added.
        prop_kw: mappable
            Can provide additional IOAPI (or other) properties for the file.
        """
        now = datetime.now()
        if hasattr(path, 'read'):
            gdf = path
            path = '<inline>'
        elif os.path.exists(path):
            gdf = open(path, 'r')
        else:
            gdf = io.StringIO(path)
            path = '<inline>'
        gdlines = gdf.read().strip().split('\n')
        gdlines = [line.strip() for line in gdlines]
        assert(gdlines[0] == "' '")
        assert(gdlines[-1] == "' '")
        i = 0
        prj = {}
        grd = OrderedDict()
        blanks = []
        while i < len(gdlines):
            line = gdlines[i]
            if line.strip() == "' '":
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
        self.GDNAM = GDNAM
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
        self.SDATE = -635
        self.STIME = 0
        self.TSTEP = 0
        self._synthvars = None
        # user supplied properties should overwrite
        # defaults above
        for pk, pv in prop_kw.items():
            self.setncattr(pk, pv)

        self.addtime()
        self.setgrid(withcf=withcf)
        self.updatemeta()
        self.getVarlist()

    def setgrid(self, key=None, withcf=True):
        if key is None:
            if self.GDNAM is None:
                self.GDNAM = key = list(self._grd)[0]
            else:
                key = self.GDNAM
        grd = self._grd[key]
        prj = self._prj[grd['PRJNAME']]
        self.setncatts(prj)
        self.setncatts(grd)
        self.adddims(withcf=withcf)

    def addtime(self):
        self.createDimension('TSTEP', 1).setunlimited(True)
        self.createDimension('DATE-TIME', 2)
        self.createDimension('VAR', 1)
        self.NVARS = 1
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
        self.createVariable(
            'DUMMY', 'f', dims,
            units='none'.ljust(16),
            long_name='DUMMY'.ljust(16),
            var_desc='DUMMY'.ljust(80)
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
