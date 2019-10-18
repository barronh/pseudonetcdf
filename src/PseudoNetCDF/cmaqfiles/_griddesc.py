import PseudoNetCDF as pnc
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
        self, path, GDNAM=None, VGLVLS=[1., 0.], VGTOP=5000., FTYPE=1, VGTYP=7
    ):
        now = datetime.now()
        gdlines = open(path, 'r').read().strip().split('\n')
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
        self.VGTOP = VGTOP
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
        self.addtime()
        self.setgrid()

    def setgrid(self, key=None):
        if key is None:
            if self.GDNAM is None:
                self.GDNAM = key = list(self._grd)[0]
            else:
                key = self.GDNAM
        grd = self._grd[key]
        prj = self._prj[grd['PRJNAME']]
        self.setncatts(prj)
        self.setncatts(grd)
        self.adddims()

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
        tflag[:, :, 0] = -635
        tflag[:, :, 1] = 0

    def adddims(self):
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
            units='none',
            long_name='DUMMY'.ljust(16),
            var_desc='DUMMY'.ljust(80)
        )
        setattr(self, 'VAR-LIST', 'DUMMY'.ljust(16))

        if self._synthvars is None:
            oldkeys = set(self.variables)
        else:
            for key in self._synthvars:
                if key in self.variables:
                    del self.variables[key]
        pnc.conventions.ioapi.add_cf_from_ioapi(self)
        if self._synthvars is None:
            self._synthvars = set(self.variables).difference(oldkeys)


if __name__ == '__main__':
    f = griddesc('GRIDDESC', '36US3')
