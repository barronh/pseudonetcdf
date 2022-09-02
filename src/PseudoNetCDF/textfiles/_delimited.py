from __future__ import print_function  # , unicode_literals
from PseudoNetCDF.sci_var import PseudoNetCDFFile
import numpy as np
from warnings import warn
from PseudoNetCDF._getwriter import registerwriter
import unittest

_coordkeys = ("time time_bounds TFLAG ETFLAG latitude latitude_bounds " +
              "longitude longitude_bounds lat lat_bnds lon lon_bnds " +
              "etam_pressure etai_pressure layer_bounds layer47 " +
              "layer").split()


class csv(PseudoNetCDFFile):
    def __init__(self, path, coordkeys=_coordkeys, delimiter=',',
                 names=True, backend=None, defaultcoord='record', **kwds):
        """
        path : str
            place to find csv file
        coordkeys : iterable of strings
            use these keys as dimensions and coordinate variables
        delimiter : str
            use this as delimiter (default = ',') renamed as sep for pandas
        names : iterable of strings or True
            with pandas backend and names==True: header='infer'
            otherwise, passed directly as keyword
        backend : str or None
            'numpy' numpy.recfromtxt or 'pandas' pandas.read_csv; defaults to
            pandas if available
        defaultcoord : str
            if no coordkeys are found, use this str to create a arbitrary
            coordinate based on each record
        kwds : mappable
            corresponds to numpy.recfromtxt or pandas.read_csv keywords

        * Note: currently only works when all coordinate variables are 1-d
        """
        if backend != 'numpy':
            try:
                import pandas
                if backend is None:
                    backend = 'pandas'
            except Exception:
                if backend == 'pandas':
                    raise ValueError(
                        'pandas library not available, try another backend'
                    )

        if backend == 'numpy':
            npkwds = kwds.copy()
            npkwds['names'] = names
            npkwds['delimiter'] = delimiter
            data = np.recfromtxt(path, **npkwds)
        elif backend == 'pandas':
            pdkwds = kwds.copy()
            if names is True:
                pdkwds.setdefault('header', 'infer')
            else:
                pdkwds['names'] = names
            pdkwds.setdefault('sep', delimiter)
            odata = pandas.read_csv(path, **pdkwds)
            # pandas leaves whitespace in names, which is not good for
            # netcdf-like names and probably not intended.
            # odata.rename(columns=lambda x: x.strip(), inplace=True)
            # relying on user to supply skipinitialspace=True
            data = odata.to_records(index=False)
        else:
            raise ValueError(
                "backend options are 'numpy' or 'pandas': got '%s'" % backend
            )

        varkeys = [vk for vk in data.dtype.names if vk not in coordkeys]
        dimkeys = tuple([dk for dk in coordkeys if dk in data.dtype.names])
        dimvars = {}
        if len(dimkeys) == 0:
            dimkeys = (defaultcoord,)
            dv = self.createDimension(defaultcoord, data.shape[0])
            dvar = self.createVariable(defaultcoord, 'd', (defaultcoord,))
            dvar[:] = dimvars[defaultcoord] = np.arange(len(dv))
        else:
            for dk in dimkeys:
                dimvars[dk] = data[dk]
                dv = np.unique(data[dk])
                dv.sort()
                self.createDimension(dk, len(dv))
                mydtype = dv.dtype.char
                if mydtype == 'S':
                    mydtype = dv.dtype
                dvar = self.createVariable(dk, mydtype, (dk,))
                dvar[:] = dv

        for vk in varkeys:
            vv = data[vk]
            if vv.dtype.char != 'S':
                var = self.createVariable(
                    vk, vv.dtype.char, dimkeys, fill_value=-999)
                var[:] = -999
            else:
                var = self.createVariable(vk, vv.dtype.char, dimkeys)

        bigidx = []
        bigidx = dict([(dk, (dimvars[dk][:, None] == self.variables[dk]
                             [None, :]).argmax(1)) for dk in dimkeys])
        for vk in varkeys:
            vv = data[vk]
            ov = self.variables[vk]
            myidx = tuple([bigidx[dk] for dk in ov.dimensions])
            ov[myidx] = vv


def ncf2csv(ifile, outpath, delimiter=',', coordkeys=_coordkeys):
    header = [k for k, v in ifile.variables.items(
    ) if k not in coordkeys and v.size > 1 and k not in ifile.dimensions]
    dims = set([ifile.variables[k].dimensions for k in header])
    if len(dims) > 1:
        if hasattr(outpath, 'write'):
            warn('Multiple csv outputs will be separated by ### because ' +
                 'not all output variables have the same dimensions')
        else:
            warn('Making multiple csv outputs because not all output ' +
                 'variables have the same dimensions')
    dimsets = {}
    for di, dim in enumerate(dims):
        if len(dims) > 1:
            if not hasattr(outpath, 'write'):
                outfile = open(outpath + str(di), 'wt')
            else:
                print('###', file=outfile)
        else:
            if not hasattr(outpath, 'write'):
                outfile = open(outpath, mode='wt')
            else:
                outfile = outpath

        dimheader = [k for k in dim if k in ifile.variables]
        header = dimsets[dim] = [
            k for k, v in ifile.variables.items() if v.dimensions == dim]
        dimvars = [ifile.variables[k] for k in dimheader]
        vars = [ifile.variables[k] for k in header]
        outtext = delimiter.join(dimheader + header)
        print(outtext, file=outfile)
        for idx in np.ndindex(ifile.variables[header[-1]].shape):
            outvals = []
            for dk, dv in zip(dimheader, dimvars):
                dv = ifile.variables[dk]
                didx = tuple([iidx for i, iidx in enumerate(
                    idx) if dim[i] in dv.dimensions])
                outvals.append(repr(dv[didx]))
            for vk, vv in zip(header, vars):
                outvals.append(repr(vv[idx]))
            outtext = delimiter.join(outvals)
            print(outtext, file=outfile)


registerwriter('csv', ncf2csv)


class TestCsv(unittest.TestCase):
    def setUp(self):
        from PseudoNetCDF import PseudoNetCDFFile
        self.checkval = """time,layer,latitude,longitude,test
0.0,0.0,0.0,0.0,0.0
0.0,0.0,0.0,1.0,1.0
0.0,0.0,0.0,2.0,2.0
0.0,0.0,0.0,3.0,3.0
0.0,0.0,0.0,4.0,4.0
0.0,0.0,1.0,0.0,5.0
0.0,0.0,1.0,1.0,6.0
0.0,0.0,1.0,2.0,7.0
0.0,0.0,1.0,3.0,8.0
0.0,0.0,1.0,4.0,9.0
0.0,0.0,2.0,0.0,10.0
0.0,0.0,2.0,1.0,11.0
0.0,0.0,2.0,2.0,12.0
0.0,0.0,2.0,3.0,13.0
0.0,0.0,2.0,4.0,14.0
0.0,0.0,3.0,0.0,15.0
0.0,0.0,3.0,1.0,16.0
0.0,0.0,3.0,2.0,17.0
0.0,0.0,3.0,3.0,18.0
0.0,0.0,3.0,4.0,19.0
0.0,1.0,0.0,0.0,20.0
0.0,1.0,0.0,1.0,21.0
0.0,1.0,0.0,2.0,22.0
0.0,1.0,0.0,3.0,23.0
0.0,1.0,0.0,4.0,24.0
0.0,1.0,1.0,0.0,25.0
0.0,1.0,1.0,1.0,26.0
0.0,1.0,1.0,2.0,27.0
0.0,1.0,1.0,3.0,28.0
0.0,1.0,1.0,4.0,29.0
0.0,1.0,2.0,0.0,30.0
0.0,1.0,2.0,1.0,31.0
0.0,1.0,2.0,2.0,32.0
0.0,1.0,2.0,3.0,33.0
0.0,1.0,2.0,4.0,34.0
0.0,1.0,3.0,0.0,35.0
0.0,1.0,3.0,1.0,36.0
0.0,1.0,3.0,2.0,37.0
0.0,1.0,3.0,3.0,38.0
0.0,1.0,3.0,4.0,39.0
0.0,2.0,0.0,0.0,40.0
0.0,2.0,0.0,1.0,41.0
0.0,2.0,0.0,2.0,42.0
0.0,2.0,0.0,3.0,43.0
0.0,2.0,0.0,4.0,44.0
0.0,2.0,1.0,0.0,45.0
0.0,2.0,1.0,1.0,46.0
0.0,2.0,1.0,2.0,47.0
0.0,2.0,1.0,3.0,48.0
0.0,2.0,1.0,4.0,49.0
0.0,2.0,2.0,0.0,50.0
0.0,2.0,2.0,1.0,51.0
0.0,2.0,2.0,2.0,52.0
0.0,2.0,2.0,3.0,53.0
0.0,2.0,2.0,4.0,54.0
0.0,2.0,3.0,0.0,55.0
0.0,2.0,3.0,1.0,56.0
0.0,2.0,3.0,2.0,57.0
0.0,2.0,3.0,3.0,58.0
0.0,2.0,3.0,4.0,59.0
1.0,0.0,0.0,0.0,60.0
1.0,0.0,0.0,1.0,61.0
1.0,0.0,0.0,2.0,62.0
1.0,0.0,0.0,3.0,63.0
1.0,0.0,0.0,4.0,64.0
1.0,0.0,1.0,0.0,65.0
1.0,0.0,1.0,1.0,66.0
1.0,0.0,1.0,2.0,67.0
1.0,0.0,1.0,3.0,68.0
1.0,0.0,1.0,4.0,69.0
1.0,0.0,2.0,0.0,70.0
1.0,0.0,2.0,1.0,71.0
1.0,0.0,2.0,2.0,72.0
1.0,0.0,2.0,3.0,73.0
1.0,0.0,2.0,4.0,74.0
1.0,0.0,3.0,0.0,75.0
1.0,0.0,3.0,1.0,76.0
1.0,0.0,3.0,2.0,77.0
1.0,0.0,3.0,3.0,78.0
1.0,0.0,3.0,4.0,79.0
1.0,1.0,0.0,0.0,80.0
1.0,1.0,0.0,1.0,81.0
1.0,1.0,0.0,2.0,82.0
1.0,1.0,0.0,3.0,83.0
1.0,1.0,0.0,4.0,84.0
1.0,1.0,1.0,0.0,85.0
1.0,1.0,1.0,1.0,86.0
1.0,1.0,1.0,2.0,87.0
1.0,1.0,1.0,3.0,88.0
1.0,1.0,1.0,4.0,89.0
1.0,1.0,2.0,0.0,90.0
1.0,1.0,2.0,1.0,91.0
1.0,1.0,2.0,2.0,92.0
1.0,1.0,2.0,3.0,93.0
1.0,1.0,2.0,4.0,94.0
1.0,1.0,3.0,0.0,95.0
1.0,1.0,3.0,1.0,96.0
1.0,1.0,3.0,2.0,97.0
1.0,1.0,3.0,3.0,98.0
1.0,1.0,3.0,4.0,99.0
1.0,2.0,0.0,0.0,100.0
1.0,2.0,0.0,1.0,101.0
1.0,2.0,0.0,2.0,102.0
1.0,2.0,0.0,3.0,103.0
1.0,2.0,0.0,4.0,104.0
1.0,2.0,1.0,0.0,105.0
1.0,2.0,1.0,1.0,106.0
1.0,2.0,1.0,2.0,107.0
1.0,2.0,1.0,3.0,108.0
1.0,2.0,1.0,4.0,109.0
1.0,2.0,2.0,0.0,110.0
1.0,2.0,2.0,1.0,111.0
1.0,2.0,2.0,2.0,112.0
1.0,2.0,2.0,3.0,113.0
1.0,2.0,2.0,4.0,114.0
1.0,2.0,3.0,0.0,115.0
1.0,2.0,3.0,1.0,116.0
1.0,2.0,3.0,2.0,117.0
1.0,2.0,3.0,3.0,118.0
1.0,2.0,3.0,4.0,119.0
"""
        testfile = self.testfile = PseudoNetCDFFile()
        testfile.createDimension('time', 2)
        testfile.createDimension('layer', 3)
        testfile.createDimension('latitude', 4)
        testfile.createDimension('longitude', 5)
        for dk, dv in testfile.dimensions.items():
            var = testfile.createVariable(dk, 'f', (dk,))
            var[:] = np.arange(len(dv), dtype='f')
        var = testfile.createVariable(
            'test', 'f', ('time', 'layer', 'latitude', 'longitude'))
        var[:] = np.arange(2 * 3 * 4 * 5).reshape(2, 3, 4, 5)

    def testNCF2CSV(self):
        from PseudoNetCDF.pncgen import pncgen
        from io import BytesIO
        import tempfile
        out = BytesIO()
        out = tempfile.TemporaryFile(mode='w+t')
        pncgen(self.testfile, out, inmode='r',
               outmode='w', format='csv', verbose=0)
        out.seek(0, 0)
        outval = out.read()
        assert (outval == self.checkval)
