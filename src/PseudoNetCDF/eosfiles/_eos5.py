__all__ = ['eos5']

from collections import OrderedDict
from ..core._files import PseudoNetCDFFile
from ..core._variables import PseudoNetCDFVariable


class _dummyvar(object):
    def __init__(self, key, var, dimensions):
        self._var = var
        for k in var.ncattrs():
            setattr(self, k, getattr(var, k))
        self.dimensions = dimensions
        self.dtype = var.dtype
        self.shape = var.shape
        self.name = key
        self.typecode = var.dtype.char

    def array(self):
        return self._var[:]

    def __repr__(self):
        return PseudoNetCDFVariable.__repr__(self)

    def __str__(self):
        return PseudoNetCDFVariable.__str__(self)

    def __getitem__(self, *args, **kwds):
        return self._var.__getitem__(*args, **kwds)

    def ncattrs(self):
        return self._var.ncattrs()

    def getncatts(self):
        return {k: getattr(self._var, k) for k in self._var.ncattrs()}


def _keyval(s):
    i = s.find('=')
    if i == -1:
        key = s
        val = ''
    else:
        key = s[:i].lstrip()
        val = s[i + 1:]
        try:
            val = eval(val)
        except Exception:
            pass

    return key, val


class eos5(PseudoNetCDFFile):
    def __init__(self, path, *args, **kwds):
        """
        Fixes unnamed dimensions for

        Arguments
        ---------
        path : str
            path to hdf-eos5 file
        args : iterable
            arguments for pncopen
        kwds : mappable
            keyword arguments for pncopen
        """
        from .. import pncopen
        self._tmpf = pncopen(path, *args, format='netcdf', **kwds)
        self.groups['HDFEOS'] = PseudoNetCDFFile()
        self.groups['HDFEOS'].groups['SWATHS'] = PseudoNetCDFFile()
        self.parse_struct()

    def parse_struct(self):
        structv = self._tmpf['HDFEOS INFORMATION/StructMetadata.0'][:]
        lines = structv.split('\n')
        linei = iter(lines)
        for line in linei:
            key, val = _keyval(line)
            if key == 'GROUP' and val.startswith('SWATH_'):
                swathf = PseudoNetCDFFile()
                swathf.groups['Geolocation Fields'] = PseudoNetCDFFile()
                swathf.groups['Data Fields'] = PseudoNetCDFFile()
                # tmp1f = PseudoNetCDFFile()
                useddims = set()
            elif key == 'SwathName':
                SwathName = swathf.SwathName = val
                self['HDFEOS/SWATHS'].groups[SwathName] = swathf
            elif key == 'DimensionName':
                dkey = val
                dummy, val = _keyval(next(linei))
                assert(dummy == 'Size')
                swathf.createDimension(dkey, int(val))
                # tmp1f.createDimension(dkey, 1)
            elif (
                key == 'OBJECT' and (
                    val.startswith('GeoField_') or
                    val.startswith('DataField_')
                )
            ):
                fldkey, fldname = _keyval(next(linei))
                _keyval(next(linei))
                dimkey, dims = _keyval(next(linei))
                if isinstance(dims, str):
                    dims = (dims,)
                _keyval(next(linei))
                startkey = 'HDFEOS/SWATHS/' + SwathName
                if val.startswith('Geo'):
                    grpkey = startkey + '/Geolocation Fields'
                elif val.startswith('Data'):
                    grpkey = startkey + '/Data Fields'
                tmpvar = self._tmpf[grpkey].variables[fldname]
                # var = tmp1f.createVariable(fldname, tmpvar.dtype, dims)
                # var.setncatts(
                #     dict([
                #         (k, tmpvar.getncattr(k))
                #         for k in tmpvar.ncattrs()
                #     ])
                # )
                tmpvar = self[grpkey].variables[fldname] = _dummyvar(fldname, tmpvar, dims)
                for di, dimkey in enumerate(dims):
                    dim = swathf.dimensions[dimkey]
                    useddims.add(dimkey)
                    if dim.isunlimited():
                        size = tmpvar.shape[di]
                        newdim = swathf.createDimension(dimkey, size)
                        newdim.isunlimited(True)
            elif key == 'END_GROUP' and val.startswith('SWATH_'):
                for dk in list(swathf.dimensions):
                    if dk not in useddims:
                        del swathf.dimensions[dk]
                        # del tmp1f.dimensions[dk]
                self.groups['HDFEOS'].groups[SwathName] = swathf


if __name__ == '__main__':
    f = eos5(
        '/work/ROMO/users/bhenders/obs/OMNO2d/aura.gesdisc.eosdis.nasa.gov/' +
        'data/Aura_OMI_Level2/OMNO2.003/2016/001/' +
        'OMI-Aura_L2-OMNO2_2016m0101t0059-o60970_v003-2018m0617t014243.he5'
    )
