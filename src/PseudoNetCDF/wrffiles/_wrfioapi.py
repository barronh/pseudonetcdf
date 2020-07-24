__all__ = ['wrf_base', 'wrf']
from ..core._files import PseudoNetCDFFile, netcdf


class wrf_base(PseudoNetCDFFile):
    @classmethod
    def isMine(self, path, *args, **kwds):
        return False

    def getTimes(self):
        cf = self.getCoordFile()
        return PseudoNetCDFFile.getTimes(cf)

    def getCoordFile(self):
        from ..conventions.ioapi import add_cf_from_wrfioapi
        cf = self.subsetVariables(list(self.getCoords()))
        add_cf_from_wrfioapi(cf)
        coordkeys = [k for k in cf.variables if k in cf.dimensions]
        coordkeys.extend(['longitude', 'latitude', 'time'])
        cf.setCoords(coordkeys)
        return cf

    def __init__(self, *args, **kwds):
        super().__init__(self, *args, **kwds)
        self._defCoords()

    def _defCoords(self):
        coordkeys = [
            key for key in self.variables
            if (
                key.startswith('XLAT') or key.startswith('XLONG') or
                key.startswith('Times')
            )
        ]
        self.setCoords(coordkeys)

    def plot(self, varkey, **kwds):
        var = self.variables[varkey]
        mykwds = kwds.copy()
        map_kw = kwds.pop('map_kw', None)
        mykwds['map_kw'] = False

        dimopts = [
            dk for dk in var.dimensions
            if len(self.dimensions[dk]) > 0
        ]
        plottype = '-'.join(dimopts[-2:][::-1])
        mykwds.setdefault('plottype', plottype)

        cf = self.getCoordFile()
        cf.copyVariable(var, key=varkey)
        ax = PseudoNetCDFFile.plot(cf, varkey, **kwds)
        maptypes = ('west_east-south_north', 'west_east_stag-south_north_stag')

        if plottype in maptypes and map_kw is not False:
            if map_kw is True or map_kw is None:
                map_kw = {}
            mymap_kw = map_kw.copy()
            mymap_kw.pop('coastlines', True)
            mymap_kw.pop('countries', True)
            mymap_kw.pop('states', False)
            mymap_kw.pop('counties', False)
            bmap = self.getMap(**mymap_kw)
            if map_kw.get('coastlines', True):
                bmap.drawcoastlines(ax=ax)
            if map_kw.get('countries', True):
                bmap.drawcountries(ax=ax)
            if map_kw.get('states', False):
                bmap.drawstates(ax=ax)
            if map_kw.get('counties', False):
                bmap.drawcounties(ax=ax)

        return ax

    def val2idx(self, *args, **kwds):
        cf = self.getCoordFile()
        return PseudoNetCDFFile.val2idx(cf, *args, **kwds)

    def xy2ll(self, *args, **kwds):
        cf = self.getCoordFile()
        return PseudoNetCDFFile.xy2ll(cf, *args, **kwds)

    def ll2xy(self, *args, **kwds):
        cf = self.getCoordFile()
        return PseudoNetCDFFile.ll2xy(cf, *args, **kwds)

    def ll2ij(self, *args, **kwds):
        cf = self.getCoordFile()
        return PseudoNetCDFFile.ll2ij(cf, *args, **kwds)

    def ij2ll(self, *args, **kwds):
        cf = self.getCoordFile()
        return PseudoNetCDFFile.ij2ll(cf, *args, **kwds)

    def getMap(self, maptype='basemap_auto', **kwds):
        mykwds = kwds.copy()
        cf = self.getCoordFile()
        if maptype.endswith('_auto'):
            x = cf.variables['x']
            y = cf.variables['y']
            xmax = x.max() - x.min()
            ymax = y.max() - y.min()
            lllon, lllat = cf.xy2ll(0, 0)
            urlon, urlat = cf.xy2ll(xmax, ymax)
            mykwds.setdefault('llcrnrlon', lllon)
            mykwds.setdefault('llcrnrlat', lllat)
            mykwds.setdefault('urcrnrlon', urlon)
            mykwds.setdefault('urcrnrlat', urlat)
            maptype = maptype[:-5]

        return PseudoNetCDFFile.getMap(cf, maptype=maptype, **mykwds)


class wrf(wrf_base, netcdf):
    def __init__(self, *args, **kwds):
        netcdf.__init__(self, *args, **kwds)
        self._defCoords()

    def _newlike(self):
        if self.get_dest() is not None:
            outf = wrf(**self.get_dest())
        elif isinstance(self, PseudoNetCDFFile):
            outt = wrf_base
            outf = outt.__new__(outt)
        else:
            outf = PseudoNetCDFFile()
        outf.set_varopt(**self.get_varopt())
        return outf

    @classmethod
    def from_ncf(cls, infile):
        outf = wrf_base()
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
            dimkeys = ['Time', 'bottom_top', 'west_east', 'north_south']
            for dimkey in dimkeys:
                if dimkey not in f.dimensions:
                    return False
            else:
                return True
        except Exception:
            return False
