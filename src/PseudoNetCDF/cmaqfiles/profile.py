from __future__ import print_function
__all__ = ['bcon_profile', 'icon_profile']
from PseudoNetCDF import PseudoNetCDFFile, PseudoNetCDFVariables
from ._ioapi import ioapi_base
from matplotlib.mlab import csv2rec
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
from datetime import datetime
import numpy as np

def _getunit(varkey):
    if varkey.startswith('A') and varkey[-1:] in ('I', 'J', 'K'):
        unit = 'micrograms/m**3'
    elif varkey.startswith('ASEA'):
        unit = 'micrograms/m**3'
    elif varkey.startswith('NUM'):
        unit = '#/m**3'
    elif varkey.startswith('SRF'):
        unit =  'm**2/m**3'
    else:
        unit = 'ppmV'
    return unit.ljust(16)
    
class icon_profile(ioapi_base):
    def __init__(self, path):
        lines = open(path).read().split('\n')
        header = lines[3].split()
        nlay, nspc = [int(_v) for _v in header[:2]]
        sigmas = [float(_v) for _v in header[2:]]
        nsigmas = len(sigmas)
        try:
            dateo = datetime.strptime(lines[4].strip(), '%Y-%m-%d')
        except:
            date, time = [int(_v) for _v in lines[4].split()]
        starts =  [5]
        ends = [s + nspc for s in starts]
        keys = ['all']
        fieldnames = ('name',) + tuple(['s%f' % i for i in sigmas])
        data = dict([(k, csv2rec(StringIO('\n'.join(lines[s:e])), delimiter = ' ', names = fieldnames, converterd = dict(names = lambda x: str(x).strip()))) for k, s, e in zip(keys, starts, ends)])
        profile_spcs = np.char.strip(data[keys[0]]['name'])
        data_type = data[keys[0]].dtype
        data_shape =  data[keys[0]].shape
        self.createDimension('sigma', nsigmas)
        self.createDimension('LAY', nsigmas - 1)
        self.createVariable('sigma', 'f', ('sigma',), values = np.array(sigmas), units = 'sigma')
        self.createVariable('LAY', 'f', ('LAY',), values = np.array(sigmas).repeat(2, 0)[1:-1].reshape(-1, 2).mean(1), units = 'sigma')
        self.VGLVLS = self.variables['sigma']
        self.VGTOP = 5000
        ks = keys[1:]
        for k in ks:
            try:
                assert((np.char.strip(data[k]['name']) == profile_spcs).all())
                assert(data[k].dtype == data_type)
                assert(data[k].dtype == data_type)
            except AssertionError:
                raise IOError('File is corrupt or inconsistent')
        varlist = []
        for a in data['all']:
            varkey = a[0].strip()
            self.createVariable(varkey, 'f', ('LAY',), units = _getunit(varkey), values = np.array([tuple(a)[1:]])[0].astype('f'), long_name = varkey.ljust(16), var_desc = varkey.ljust(16))
            varlist.append(varkey.ljust(16))
        self.NVARS = len(varlist)
        self.createDimension('VAR', self.NVARS)
        setattr(self, 'VAR-LIST', ''.join(varlist))
            
    
    def toGrid(self, ncols, nrows, nsteps = 1):
        """
        Parameters
        ----------
        ncols  : number of ncols for output file
        nrows  : number of rows for output file
        nsteps : number of output time steps
        """
        out = self.subsetVariables(['sigma', 'LAY'], exclude = True)
        out = out.insertDimension(TSTEP = nsteps, before = 'LAY')
        out = out.insertDimension(ROW = nrows, after = 'LAY')
        out = out.insertDimension(COL = ncols, after = 'ROW')
        return out
        
    def toGriddedFile(self, ioapifile):
        """
        Parameters
        ----------
        ioapifile : file with IOAPI metadata
        """
        nsteps = len(ioapifile.dimensions['TSTEP'])
        nrows = len(ioapifile.dimensions['ROW'])
        ncols = len(ioapifile.dimensions['COL'])
        # interp to ioapi vglvls
        out = self.interpSigma(vglvls = ioapifile.VGLVLS, vgtop = ioapifile.VGTOP)
        # repeat to grid
        out = out.toGrid(ncols, nrows, nsteps)
        for pk in ioapifile.ncattrs():
            if pk in ('VAR-LIST', 'FILEDESC', 'HISTORY'): continue
            setattr(out, pk , getattr(ioapifile, pk))
        out.updatemeta()
        return out
        
class bcon_profile(ioapi_base):
    def __init__(self, path):
        lines = open(path).read().split('\n')
        header = lines[3].split()
        nlay, nspc = [int(_v) for _v in header[:2]]
        sigmas = [float(_v) for _v in header[2:]]
        nsigmas = len(sigmas)
        try:
            dateo = datetime.strptime(lines[4].strip(), '%Y-%m-%d')
        except:
            date, time = [int(_v) for _v in lines[4].split()]
        starts =  [5 + i + i * nspc for i in range(4)]
        ends = [s + 1 + nspc for s in starts]
        keys = [lines[s].strip().lower() for s in starts]
        fieldnames = ('name',) + tuple(['s%f' % i for i in sigmas])
        data = dict([(k, csv2rec(StringIO('\n'.join(lines[s+1:e])), delimiter = ' ', names = fieldnames, converterd = dict(names = lambda x: str(x).strip()))) for k, s, e in zip(keys, starts, ends)])
        profile_spcs = np.char.strip(data[keys[0]]['name'])
        data_type = data[keys[0]].dtype
        data_shape =  data[keys[0]].shape
        self.createDimension('sigma', nsigmas)
        self.createDimension('LAY', nsigmas - 1)
        self.createDimension('south_east_north_west', 4)
        self.createVariable('sigma', 'f', ('sigma',), values = np.array(sigmas), units = 'sigma')
        self.createVariable('LAY', 'f', ('LAY',), values = np.array(sigmas).repeat(2, 0)[1:-1].reshape(-1, 2).mean(1), units = 'sigma')
        ks = keys[1:]
        self.VGLVLS = self.variables['sigma']
        self.VGTOP = 5000
        for k in ks:
            try:
                assert((np.char.strip(data[k]['name']) == profile_spcs).all())
                assert(data[k].dtype == data_type)
                assert(data[k].dtype == data_type)
            except AssertionError:
                raise IOError('File is corrupt or inconsistent')
        varlist = []
        for w, s, e, n in zip(data['west'], data['south'], data['east'], data['north']):
            assert(w[0] == s[0] and e[0] == n[0] and n[0] == s[0])
            varkey = w[0].strip()
            self.createVariable(varkey, 'f', ('LAY', 'south_east_north_west'), units = _getunit(varkey), values = np.array([(lambda x: tuple(x)[1:])(_v) for _v in[s, e, n, w]]).T.astype('f'), long_name = varkey.ljust(16), var_desc = varkey.ljust(16))
            varlist.append(varkey.ljust(16))
        self.NVARS = len(varlist)
        self.createDimension('VAR', self.NVARS)
        setattr(self, 'VAR-LIST', ''.join(varlist))

        
    def toPerim(self, ncols, nrows, nsteps = 1):
        outf = ioapi_base()
        nperim = nrows * 2 + ncols * 2 + 4
        outf.createDimension('TSTEP', nsteps)
        outf.createDimension('DATE-TIME', 2)
        outf.createDimension('LAY', len(self.dimensions['LAY']))
        outf.createDimension('VAR', len(self.dimensions['VAR']))
        outf.createDimension('PERIM', nperim)
        varlist = []
        for vk, vv in self.variables.items():
            if vk in ('sigma', 'LAY'): continue
            outv = outf.copyVariable(vv, key = vk, dtype = 'f', dimensions = ('TSTEP', 'LAY', 'PERIM'), withdata = False)
            s = start = 0
            e = end = ncols + 1
            outv[:, :, s:e] = vv[None, :,[0]]
            s = e
            e += nrows + 1
            outv[:, :, s:e] = vv[None, :,[1]]
            s = e
            e += ncols + 1
            outv[:, :, s:e] = vv[None, :,[2]]
            s = e
            e += nrows + 1
            outv[:, :, s:e] = vv[None, :,[3]]
            varlist.append(vk.ljust(16))
        outf.NVARS = len(varlist)
        setattr(outf, 'VAR-LIST', ''.join(varlist))
        outf.updatemeta()
        return outf
        
    def toPerimFile(self, ioapifile):
        """
        Parameters
        ----------
        ioapifile : file with IOAPI metadata
        """
        nsteps = len(ioapifile.dimensions['TSTEP'])
        try:
            nrows = len(ioapifile.dimensions['ROW'])
            ncols = len(ioapifile.dimensions['COL'])
        except:
            nrows = ioapifile.NROWS
            ncols = ioapifile.NCOLS
        # interp to ioapi vglvls
        out = self.interpSigma(vglvls = ioapifile.VGLVLS, vgtop = ioapifile.VGTOP)
        # repeat to grid
        out = out.toPerim(ncols, nrows, nsteps)
        for pk in ioapifile.ncattrs():
            if pk in ('VAR-LIST', 'FILEDESC', 'HISTORY'): continue
            setattr(out, pk , getattr(ioapifile, pk))
        out.updatemeta()
        return out        


if __name__ == '__main__':
    po = profile('testdata/profile.dat')
    from PseudoNetCDF.pncdump import pncdump
    pncdump(po)
    print(po.variables['ATOL1J'].shape)
    print(po.variables['ATOL1J'].dimensions)
