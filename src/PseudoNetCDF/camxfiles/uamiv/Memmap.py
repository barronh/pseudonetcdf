__all__ = ['uamiv']
__doc__ = """
.. _Memmap
:mod:`Memmap` -- uamiv Memmap interface
============================================

.. module:: Memmap
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` memory map for CAMx
              uamiv files.  See PseudoNetCDF.sci_var.PseudoNetCDFFile
              for interface details
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""

# Distribution packages
import unittest
# from warnings import warn

# Site-Packages
from numpy import array, memmap, dtype, linspace

# This Package modules
from PseudoNetCDF.sci_var import PseudoIOAPIVariable
from PseudoNetCDF.sci_var import PseudoNetCDFVariables
from PseudoNetCDF.ArrayTransforms import ConvertCAMxTime
from PseudoNetCDF.camxfiles.units import get_uamiv_units, get_chemparam_names
from PseudoNetCDF.cmaqfiles import ioapi_base
# from PseudoNetCDF.sci_var import PseudoNetCDFFile
# from PseudoNetCDF.conventions.ioapi import add_cf_from_ioapi
# for use in identifying uncaught nan


class uamiv(ioapi_base):
    """
    uamiv provides a PseudoNetCDF interface for CAMx
    uamiv files.  Where possible, the inteface follows
    IOAPI conventions (see www.baronams.com).

    ex:
        >>> uamiv_path = 'camx_uamiv.bin'
        >>> uamivfile = uamiv(uamiv_path)
        >>> uamivfile.variables.keys()
        ['TFLAG', 'O3', 'NO', 'NO2', ...]
        >>> tflag = uamivfile.variables['TFLAG']
        >>> tflag.dimensions
        ('TSTEP', 'VAR', 'DATE-TIME')
        >>> tflag[0,0,:]
        array([2005185,       0])
        >>> tflag[-1,0,:]
        array([2005185,  240000])
        >>> v = uamivfile.variables['O3']
        >>> v.dimensions
        ('TSTEP', 'LAY', 'ROW', 'COL')
        >>> v.shape
        (25, 28, 65, 83)
        >>> uamivfile.dimensions
        {'TSTEP': 25, 'LAY': 28, 'ROW': 65, 'COL': 83}
    """

    __ione = 1
    __idum = 0
    __rdum = 0.

    def _make_header_fmt(self, ep=None):
        if ep is None:
            ep = self.__endianprefix
        self.__emiss_hdr_fmt = dtype(dict(
            names=['SPAD', 'name', 'note', 'itzon', 'nspec', 'ibdate', 'btime',
                   'iedate', 'etime', 'EPAD'],
            formats=['i', '(10, 4)S1', '(60,4)S1', 'i', 'i', 'i', 'f', 'i',
                     'f', 'i'])).newbyteorder(ep)
        self.__grid_hdr_fmt = dtype(dict(
            names=['SPAD', 'plon', 'plat', 'iutm', 'xorg', 'yorg', 'delx',
                   'dely', 'nx', 'ny', 'nz', 'iproj', 'istag', 'tlat1',
                   'tlat2', 'rdum', 'EPAD'],
            formats=['i', 'f', 'f', 'i', 'f', 'f', 'f', 'f', 'i', 'i', 'i',
                     'i', 'i', 'f', 'f', 'f', 'i'])).newbyteorder(ep)
        self.__cell_hdr_fmt = dtype(dict(
            names=['SPAD', 'ione1', 'ione2', 'nx', 'ny', 'EPAD'],
            formats=['i', 'i', 'i', 'i', 'i', 'i'])).newbyteorder(ep)
        self.__time_hdr_fmt = dtype(dict(
            names=['SPAD', 'ibdate', 'btime', 'iedate', 'etime', 'EPAD'],
            formats=['i', 'i', 'f', 'i', 'f', 'i'])).newbyteorder(ep)
        self.__spc_fmt = dtype("(10,4)S1").newbyteorder(ep)

    def __init__(self, rf, mode='r', P_ALP=None, P_BET=None, P_GAM=None,
                 XCENT=None, YCENT=None, GDTYP=None, endian='big',
                 chemparam=None):
        """
        Parameters
        ----------
        rf : string or RecordFile
            usually a path to a CAMx formatted file, can be a FortranFileUtil
            RecordFile object.
        mode : string
            file open mode read ('r'), write ('w'), append ('a', 'r+')
        P_ALP : float
            see IOAPI GRIDDESC documentation
        P_BET : float
            see IOAPI GRIDDESC documentation
        P_GAM : float
            see IOAPI GRIDDESC documentation
        XCENT : float
            see IOAPI GRIDDESC documentation
        YCENT : float
            see IOAPI GRIDDESC documentation
        GDTYP : float
            see IOAPI GRIDDESC documentation
        endian : string
            'big' or 'little' usually only if mistaken compile
        chemparam : None or string
            used to identify gases and aerosols

        Returns
        -------
        outf : uamiv
            PseudoNetCDFFile populated from file
        """
        if chemparam is None:
            self._aerosol_names = None
        else:
            self._aerosol_names = get_chemparam_names(chemparam)['aerosol']

        self.__endianprefix = dict(big='>', little='<')[endian]
        self._make_header_fmt()
        self.__rffile = rf
        self.__mode = mode

        self.createDimension('DATE-TIME', 2)

        self.__readheader()

        # Add IOAPI metavariables
        self.NLAYS = len(self.dimensions['LAY'])
        self.NROWS = len(self.dimensions['ROW'])
        self.NCOLS = len(self.dimensions['COL'])
        self.NVARS = len(self.dimensions['VAR'])
        self.NSTEPS = len(self.dimensions['TSTEP'])
        varlist = "".join([i.ljust(16) for i in self.__var_names__])
        setattr(self, 'VAR-LIST', varlist)
        self.NAME = self.__emiss_hdr['name'][0, :, 0].copy().view('S10')[
            0].decode()
        self.NOTE = self.__emiss_hdr['note'][0, :, 0].copy().view('S60')[
            0].decode()
        self.ITZON = self.__emiss_hdr['itzon'][0]
        self.FTYPE = 1
        self.VGTYP = 2
        self.VGTOP = 10000.
        self.VGLVLS = linspace(0, 1, self.NLAYS + 1)[::-1]
        self.GDNAM = "CAMx            "
        self.UPNAM = "CAMx            "
        self.FILEDESC = "CAMx            "
        # Create variables
        self.variables = PseudoNetCDFVariables(
            self.__variables, ['TFLAG', 'ETFLAG'] + self.__var_names__)
        tflag = ConvertCAMxTime(self.__memmap__['DATE']['BDATE'],
                                self.__memmap__['DATE']['BTIME'],
                                self.NVARS)
        etflag = ConvertCAMxTime(self.__memmap__['DATE']['EDATE'],
                                 self.__memmap__['DATE']['ETIME'],
                                 self.NVARS)
        tflagv = self.createVariable('TFLAG', 'i',
                                     ('TSTEP', 'VAR', 'DATE-TIME'),
                                     values=tflag, units='DATE-TIME',
                                     long_name='TFLAG'.ljust(16),
                                     var_desc='TFLAG'.ljust(80))
        etflagv = self.createVariable('ETFLAG', 'i',
                                      ('TSTEP', 'VAR', 'DATE-TIME'),
                                      values=etflag, units='DATE-TIME',
                                      long_name='ETFLAG'.ljust(16),
                                      var_desc='Ending TFLAG'.ljust(80))

        self.SDATE, self.STIME = self.variables['TFLAG'][0, 0, :]
        self.TSTEP = etflagv[0, 0, 1] - tflagv[0, 0, 1]
        if P_ALP is not None:
            self.P_ALP = P_ALP
        if P_BET is not None:
            self.P_BET = P_BET
        if P_GAM is not None:
            self.P_GAM = P_GAM
        if XCENT is not None:
            self.XCENT = XCENT
        if YCENT is not None:
            self.YCENT = YCENT
        if GDTYP is not None:
            self.GDTYP = GDTYP
        # try:
        #     add_cf_from_ioapi(self)
        # except Exception as e:
        #     warn(repr(e))
        #     pass

    def __checkfilelen(self):
        f = open(self.__rffile, 'rb')
        f.seek(0, 2)
        flen = f.tell()
        f.close()
        return flen

    @classmethod
    def isMine(cls, path):
        self = uamiv.__new__(uamiv)
        cls._make_header_fmt(self, '>')
        offset = 0
        emiss_hdr = memmap(
            path, mode='r', dtype=self.__emiss_hdr_fmt, shape=1, offset=offset)
        name = emiss_hdr['name'][0, :, 0].copy().view('S10')[
            0].decode().strip()
        if name in ('BOUNDARY', 'PTSOURCE'):
            return False

        if (
            not emiss_hdr['SPAD'] == emiss_hdr['EPAD'] and
            emiss_hdr['SPAD'] == emiss_hdr.dtype.itemsize - 8
        ):
            return False
        offset += emiss_hdr.dtype.itemsize * emiss_hdr.size

        grid_hdr = memmap(
            path, mode='r', dtype=self.__grid_hdr_fmt, shape=1, offset=offset)
        if (
            not grid_hdr['SPAD'] == grid_hdr['EPAD'] and
            grid_hdr['SPAD'] == grid_hdr.dtype.itemsize - 8
        ):
            return False
        offset += grid_hdr.dtype.itemsize * grid_hdr.size
        cell_hdr = memmap(
            path, mode='r', dtype=self.__cell_hdr_fmt, shape=1, offset=offset)
        if (
            not cell_hdr['SPAD'] == cell_hdr['EPAD'] and
            cell_hdr['SPAD'] == cell_hdr.dtype.itemsize - 8
        ):
            return False
        return name in ('AIRQUALITY', 'EMISSIONS', 'INSTANT', 'AVERAGE')

    def __readheader(self):
        ep = self.__endianprefix
        offset = 0
        self.__emiss_hdr = memmap(
            self.__rffile, mode=self.__mode, dtype=self.__emiss_hdr_fmt,
            shape=1, offset=offset)
        nspec = self.__emiss_hdr['nspec'][0]
        offset += self.__emiss_hdr.dtype.itemsize * self.__emiss_hdr.size

        self.__grid_hdr = memmap(
            self.__rffile, mode=self.__mode, dtype=self.__grid_hdr_fmt,
            shape=1, offset=offset)

        self.XORIG = self.__grid_hdr['xorg'][0]
        self.YORIG = self.__grid_hdr['yorg'][0]
        self.XCELL = self.__grid_hdr['delx'][0]
        self.YCELL = self.__grid_hdr['dely'][0]
        self.PLON = plon = self.__grid_hdr['plon'][0]
        self.PLAT = plat = self.__grid_hdr['plat'][0]
        self.TLAT1 = tlat1 = self.__grid_hdr['tlat1'][0]
        self.TLAT2 = tlat2 = self.__grid_hdr['tlat2'][0]
        self.IUTM = iutm = self.__grid_hdr['iutm'][0]
        self.ISTAG = self.__grid_hdr['istag'][0]
        self.CPROJ = cproj = self.__grid_hdr['iproj'][0]
        if hasattr(self, 'GDTYP'):
            GDTYPE = self.GDTYP
        else:
            GDTYPE = self.GDTYP = {0: 1, 1: 5, 2: 2, 3: 6}[cproj]
        if (
            cproj == 0 or
            not all([x == 0 for x in [plon, plat, tlat1, tlat2, iutm, cproj]])
        ):
            # Map CAMx projection constants to IOAPI
            self.XCENT = plon
            self.YCENT = plat
            if GDTYPE in (1, 2):
                self.P_ALP = tlat1
                self.P_BET = tlat2
                self.P_GAM = plon
            elif GDTYPE == 5:
                self.P_ALP = iutm
                self.P_BET = 0.
                self.P_GAM = 0.
            elif GDTYPE == 6:
                self.P_ALP = {90: 1, -90: -1}[plat]
                self.P_BET = tlat1
                self.P_GAM = plon
            else:
                raise ValueError('Unknown projection')

        nx = self.__grid_hdr['nx'][0]
        ny = self.__grid_hdr['ny'][0]
        nz = max(self.__grid_hdr['nz'], array([1]))[0]

        offset += self.__grid_hdr.dtype.itemsize * self.__grid_hdr.size
        self.__cell_hdr = memmap(
            self.__rffile, mode=self.__mode, dtype=self.__cell_hdr_fmt,
            shape=1, offset=offset)

        offset += self.__cell_hdr.dtype.itemsize * self.__cell_hdr.size + 4
        self.__spc_hdr = memmap(self.__rffile, mode=self.__mode,
                                dtype=self.__spc_fmt, shape=nspec,
                                offset=offset)

        offset += self.__spc_hdr.dtype.itemsize * self.__spc_hdr.size + 4

        date_time_fmt = dtype(dict(
            names=['SPAD', 'BDATE', 'BTIME', 'EDATE', 'ETIME', 'EPAD'],
            formats=['i', 'i', 'f', 'i', 'f', 'i'])).newbyteorder(ep)
        date_time_block_size = 6
        spc_1_lay_fmt = dtype(dict(
            names=['SPAD', 'IONE', 'SPC', 'DATA', 'EPAD'],
            formats=['i', 'i', '(10,4)S1',
                     '(%d,%d)f' % (ny, nx), 'i'])).newbyteorder(ep)
        spc_1_lay_block_size = 13 + nx * ny
        spc_3d_fmt = dtype((spc_1_lay_fmt, (nz,)))

        # Get species names from spc_hdr
        var_names = [spc[:, 0].copy().view('S10')[0] for spc in self.__spc_hdr]
        var_names = [v.decode() if hasattr(v, 'decode')
                     else v for v in var_names]
        self.__var_names__ = [''.join(v).strip() for v in var_names]

        data_block_fmt = dtype(dict(
            names=['DATE'] + self.__var_names__,
            formats=[date_time_fmt] + [spc_3d_fmt] * nspec))

        data_block_size = (date_time_block_size + nspec * nz *
                           spc_1_lay_block_size)
        f = open(self.__rffile)
        f.seek(0, 2)
        size = f.tell()
        f.close()
        del f
        ntimes = float(size - offset) / 4. / data_block_size
        if int(ntimes) != ntimes:
            raise ValueError(
                ("Partial time output (%f times); partial time indicates " +
                 "incomplete file.") % ntimes)
        ntimes = int(ntimes)

        self.createDimension('LAY', nz)
        self.createDimension('COL', nx)
        self.createDimension('ROW', ny)
        tstep = self.createDimension('TSTEP', ntimes)
        tstep.setunlimited(True)
        self.createDimension('VAR', nspec)

        self.__memmap__ = memmap(self.__rffile, mode=self.__mode,
                                 dtype=data_block_fmt, offset=offset)

    def __variables(self, k):
        dimensions = ('TSTEP', 'LAY', 'ROW', 'COL')
        outvals = self.__memmap__[k]['DATA']
        units = get_uamiv_units(self.NAME, k, self._aerosol_names)

        return PseudoIOAPIVariable(self, k, 'f', dimensions, values=outvals,
                                   units=units)


class TestMemmap(unittest.TestCase):
    def runTest(self):
        pass

    def setUp(self):
        pass

    def testAvg(self):
        import PseudoNetCDF.testcase
        self.assertTrue(uamiv.isMine(
            PseudoNetCDF.testcase.camxfiles_paths['uamiv']))
        emissfile = uamiv(PseudoNetCDF.testcase.camxfiles_paths['uamiv'])
        emissfile.variables['TFLAG']
        v = emissfile.variables['NO2']
        checkv = array(
            [0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
             0.00000000e+00, 0.00000000e+00, 1.24175494e-04, 2.79196858e-04,
             1.01672206e-03, 4.36782313e-04, 0.00000000e+00, 1.54810550e-04,
             3.90250643e-04, 6.18023798e-04, 3.36963218e-04, 0.00000000e+00,
             1.85579920e-04, 1.96825975e-04, 2.16468165e-04, 2.19882189e-04],
            dtype='f').reshape(1, 1, 4, 5)
        self.assertTrue((v == checkv).all())

    def testClose(self):
        import PseudoNetCDF.testcase
        emissfile = uamiv(PseudoNetCDF.testcase.camxfiles_paths['uamiv'])
        emissfile.close()

    def testSync(self):
        import PseudoNetCDF.testcase
        emissfile = uamiv(PseudoNetCDF.testcase.camxfiles_paths['uamiv'])
        emissfile.sync()


if __name__ == '__main__':
    unittest.main()
