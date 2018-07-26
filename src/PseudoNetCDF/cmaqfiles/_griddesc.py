import PseudoNetCDF as pnc
from collections import OrderedDict

"""
' '
'LamCon_40N_97W'
  2        33.000        45.000       -97.000       -97.000        40.000
' '
'SE52BENCH'
'LamCon_40N_97W'    792000.000  -1080000.000     12000.000     12000.000 100  80   1
'36US3'
'LamCon_40N_97W'  -2952000.000  -2772000.000     36000.000     36000.000 172  148   1
' '
"""
_geop = ['GDTYP', 'P_ALP', 'P_BET', 'P_GAM', 'XCENT', 'YCENT']
_prjp = ['PRJNAME', 'XORIG', 'YORIG', 'XCELL', 'YCELL', 'NCOLS', 'NROWS', 'NTHIK']
class griddesc(pnc.PseudoNetCDFFile):
    def __init__(self, path, GDNAM = None, NLAYS = 1):
        gdlines = open(path, 'r').read().strip().split('\n')
        gdlines = [l.strip() for l in gdlines]
        assert(gdlines[0] == "' '")
        assert(gdlines[-1] == "' '")
        i = 0
        geo = {}
        grd = OrderedDict()
        blanks = []
        grdnames = []
        while len(blanks) < 3:
            l = gdlines[i]
            if l.strip() == "' '":
                blanks.append(i)
            else:
                i += 1
                parts = [eval(p) for p in gdlines[i].split()]
                key = eval(l)
                if len(blanks) == 1:
                    geo[key] = dict(zip(_geop, parts))
                elif len(blanks) == 2:
                    grd[key] = dict(zip(_prjp, parts))
                else:
                   pass
            i += 1
        self._geo = geo
        self._grd = grd 
        self.GDNAM = GDNAM
        self.NLAYS = NLAYS
        self.addtime()
        self.setgrid()

    def setgrid(self, key = None):
        if key is None:
            if self.GDNAM is None:
                key = self.grd.keys()[0]
            else:
                key = self.GDNAM
        grd = self._grd[key]
        geo = self._geo[grd['PRJNAME']]
        self.setncatts(geo)
        self.setncatts(grd)
        self.adddims()

    def addtime(self):
        self.createDimension('TSTEP', 1).setunlimited(True)
        self.createDimension('DATE-TIME', 2)
        self.createDimension('VAR', 1)
        tflag = self.createVariable('TFLAG', 'i',
            ('TSTEP', 'VAR', 'DATE-TIME'),
            units='<YYYYJJJ,HHMMSS>',
            long_name='TFLAG',
            var_desc='TFLAG'
        )
        tflag[:, :, 0] = -635
        tflag[:, :, 1] = 0

    def adddims(self):
        self.createDimension('LAY', self.NLAYS)
        self.createDimension('ROW', self.NROWS)
        self.createDimension('COL', self.NCOLS)
        pnc.conventions.ioapi.add_cf_from_ioapi(self)

if __name__ == '__main__':
    f = griddesc('GRIDDESC', '36US3')
