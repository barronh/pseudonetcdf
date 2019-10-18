import numpy as np
from PseudoNetCDF import PseudoNetCDFFile
from datetime import datetime
from .arltime import arl2time

strptime = datetime.strptime

_packhdrfmt = np.dtype('>i,>S4,>i,>i')

_ijcdt = np.dtype(dict(names='IJC',
                       formats='>i2,>i2,>f'.split(',')))


class arlconcdump(PseudoNetCDFFile):
    @classmethod
    def isMine(cls, *args, **kwds):
        kwds.pop('metaonly', None)
        try:
            cls(*args, metaonly=True, **kwds)
            return True
        except Exception:
            return False

    def __init__(self, path, metaonly=False, vector=False):
        self._infile = open(path, 'rb')
        self._vector = vector
        self._infile.seek(0, 0)
        self._rec12345()
        if not metaonly:
            self._datarec()
        lat = self.createVariable(
            'latitude', 'f', ('latitude',), units='degrees_north',
            values=(self.LLCRNR_LAT + self.DELTA_LAT / 2 +
                    np.arange(self.NLATS) * self.DELTA_LAT)
        )
        lon = self.createVariable(
            'longitude', 'f', ('longitude',), units='degrees_east',
            values=(self.LLCRNR_LON + self.DELTA_LON / 2 +
                    np.arange(self.NLONS) * self.DELTA_LON)
        )
        self.createVariable(
            'latitude_bounds', 'f', ('latitude', 'nv'), units='degrees_north',
            values=np.array([lat - self.DELTA_LAT / 2,
                             lat + self.DELTA_LAT / 2]).T
        )
        self.createVariable(
            'longitude_bounds', 'f', ('longitude', 'nv'), units='degrees_east',
            values=np.array([lon - self.DELTA_LON / 2,
                             lon + self.DELTA_LON / 2]).T
        )
        if self._vector:
            tempf = self.subsetVariables(['latitude', 'latitude_bounds',
                                          'longitude', 'longitude_bounds',
                                          'time', 'time_bounds',
                                          'layer'])\
                        .sliceDimensions(latitude=self.variables['J'],
                                         longitude=self.variables['I'],
                                         time=self.variables['T'],
                                         layer=self.variables['K'])\
                        .renameDimensions(latitude='points',
                                          longitude='points',
                                          time='points',
                                          layer='points')
            for key in list(tempf.variables):
                self.variables[key] = tempf.variables[key]

    def _rec12345(self):
        """
        Record #1

        CHAR*4 Meteorological MODEL Identification
        INT*4 Meteorological file starting time (YEAR, MONTH, DAY, HOUR,
                                                 FORECAST-HOUR)
        INT*4 NUMBER of starting locations
        INT*4 Concentration packing flag (0=no 1=yes)
        """
        infile = self._infile
        rec1 = np.fromfile(
            infile, dtype='>i,>S4,>i,>i,>i,>i,>i,>i,>i,>i', count=1)[0]
        assert(rec1['f0'] == rec1['f9'])
        nloc = self.NSTARTLOCS = rec1['f7']
        rec1['f8'] == 1
        self.METMODEL = rec1['f1'].decode()
        self.MET_YEAR = rec1['f2']
        self.MET_MONTH = rec1['f3']
        self.MET_DAY = rec1['f4']
        self.MET_HOUR = rec1['f5']
        self.MET_FORECAST_HOUR = rec1['f6']
        self.STARTING_LOCATIONS = rec1['f7']
        self.PACKED = rec1['f8']

        """
        Record #2 Loop to record: Number of starting locations

        INT*4 Release starting time (YEAR, MONTH, DAY, HOUR)
        REAL*4 Starting location and height (LATITUDE, LONGITUDE, METERS)
        INT*4 Release starting time (MINUTES)
        """
        rec2 = np.fromfile(infile, dtype='>i,>4i,>3f,>i,>i', count=nloc)
        assert((rec2['f0'] == rec2['f4']).all())
        self.createDimension('starts', nloc)
        sy = self.createVariable('START_YEAR', 'i', ('starts',), units='year')
        sy[:] = rec2['f1'][:, 0]
        sm = self.createVariable(
            'START_MONTH', 'i', ('starts',), units='month')
        sm[:] = rec2['f1'][:, 1]
        sd = self.createVariable('START_DAY', 'i', ('starts',), units='day')
        sd[:] = rec2['f1'][:, 2]
        sh = self.createVariable('START_HOUR', 'i', ('starts',), units='hour')
        sh[:] = rec2['f1'][:, 3]
        slat = self.createVariable(
            'START_LAT', 'f', ('starts',), units='degrees_north')
        slat[:] = rec2['f2'][:, 0]
        slon = self.createVariable(
            'START_LON', 'f', ('starts',), units='degrees_east')
        slon[:] = rec2['f2'][:, 1]
        salt = self.createVariable(
            'START_ALT', 'f', ('starts',), units='meters')
        salt[:] = rec2['f2'][:, 2]

        """
        Record #3

        INT*4 Number of (LATITUDE-POINTS, LONGITUDE-POINTS)
        REAL*4 Grid spacing (DELTA-LATITUDE,DELTA-LONGITUDE)
        REAL*4 Grid lower left corner (LATITUDE, LONGITUDE)
        """
        rec3 = np.fromfile(infile, dtype='>i,>2i,>2f,>2f,>i', count=1)[0]
        assert(rec3['f0'] == rec3['f4'])

        self.NLATS = rec3['f1'][0]
        self.NLONS = rec3['f1'][1]
        self.DELTA_LAT = rec3['f2'][0]
        self.DELTA_LON = rec3['f2'][1]
        self.LLCRNR_LAT = rec3['f3'][0]
        self.LLCRNR_LON = rec3['f3'][1]
        """
        Record #4

        INT*4 NUMBER of vertical levels in concentration grid
        INT*4 HEIGHT of each level (meters above ground)
        """
        tmp = np.fromfile(infile, dtype='>i,>i', count=1)[0]
        nlays = self.NLAYS = tmp['f1']
        infile.seek(-8, 1)
        rec4 = np.fromfile(
            infile, dtype='>i,>i,>{}i,>i'.format(nlays), count=1)[0]
        assert(rec4['f0'] == rec4['f3'])
        self.createDimension('layer', nlays)
        var = self.createVariable('layer', 'i', ('layer',))
        var.units = 'meters agl'
        var[:] = rec4['f2']

        """
        Record #5

        INT*4 NUMBER of different pollutants in grid
        CHAR*4 Identification STRING for each pollutant
        """
        tmp = np.fromfile(infile, dtype='>i,>i', count=1)
        npols = self.NPOLS = tmp['f1'][0]
        if npols > 1 and self._vector:
            raise IOError('The vector option and multipollutant files ' +
                          'are not compatible')
        infile.seek(-8, 1)
        rec5 = np.fromfile(infile, dtype='>i,>i,>({},)S4,>i'.format(
            npols), count=1).squeeze()
        assert(rec5['f0'] == rec5['f3'])
        self.POLLUTANTS = b','.join(rec5['f2']).decode()

    def _datarec(self):
        infile = self._infile
        nlats = self.NLATS
        nlons = self.NLONS
        nlays = self.NLAYS
        npols = self.NPOLS
        pack = self.PACKED == 1
        now = infile.tell()
        infile.seek(0, 2)
        total = infile.tell()
        infile.seek(now, 0)
        datasize = total - now
        datasize
        thdr = np.dtype('>i,>6i,>i')
        nopackdfmt = '>({},{})f'.format(nlats, nlons)
        nopackfmt = np.dtype(
            dict(
                names=['B', 'POL', 'LAY', 'data', 'E'],
                formats=['>i', '>S4', '>i', nopackdfmt, '>i']
            )
        )
        outs = []
        pols = []
        lays = []
        starts = []
        stops = []
        ti = 0
        while infile.tell() != total:
            """
            Record #6 Loop to record: Number of output times

            INT*4 Sample start (YEAR MONTH DAY HOUR MINUTE FORECAST)
            """
            start = np.fromfile(infile, dtype=thdr, count=1)
            starts.append(start['f1'][0])
            """
            Record #7 Loop to record: Number of output times

            INT*4 Sample stop (YEAR MONTH DAY HOUR MINUTE FORECAST)
            """
            stop = np.fromfile(infile, dtype=thdr, count=1)
            stops.append(stop['f1'][0])
            """
            Record #8 Loop to record: Number levels, Number of pollutant types

            CHAR*4 Pollutant type identification STRING
            INT*4 Output LEVEL (meters) of this record

            No Packing (all elements)
            REAL*4 Concentration output ARRAY

            Packing (only non-zero elements)

            INT*4 Loop non-zero elements
            INT*2 First (I) index value
            INT*2 - Second (J) index value
            REAL*4 - Concentration at (I,J)**_
            """
            if pack:
                # c = np.zeros((npols, nlays, nlats, nlons), dtype='f')
                for pi in range(npols):
                    for li in range(nlays):
                        tmp = np.fromfile(infile, _packhdrfmt, count=1)
                        myp = tmp['f1'][0]
                        myl = tmp['f2'][0]
                        pols.append(myp)
                        lays.append(myl)
                        myc = tmp['f3'][0]
                        tmp = np.fromfile(infile,
                                          np.dtype([('data', _ijcdt, myc),
                                                    ('f1', '>i', 1)]), count=1)
                        tdata = tmp['data'][0]
                        Jidx = tdata['J'] - 1
                        Iidx = tdata['I'] - 1

                        # c[pi, li, Jidx, Iidx] = tdata['C']
                        if self._vector:
                            if myc > 0:
                                zl = np.array(np.ones_like(Jidx), ndmin=1)
                                outs.append([
                                    ti * zl, pi * zl, li * zl, Jidx,
                                    Iidx, tdata['C']]
                                )
                        else:
                            outs.append([ti, pi, li, Jidx, Iidx, tdata['C']])
            else:
                npblock = np.fromfile(
                    infile, dtype=nopackfmt, count=npols * nlays)
                c = npblock['data']
                outs.append(c.reshape(npols, nlays, nlats, nlons))
            ti += 1

        ntimes = ti
        if isinstance(outs[0], list):
            if self._vector:
                datablock = np.concatenate(
                    [
                        np.array(out, ndmin=2, dtype='f').reshape(6, -1)
                        for out in outs
                    ],
                    axis=1
                ).T
                npoints = datablock.shape[0]
            else:
                datablock = np.zeros((ntimes, npols, nlays, nlats, nlons),
                                     dtype='f')
                for ti, pi, li, ji, ii, v in outs:
                    datablock[ti, pi, li, ji, ii] = v
        else:
            datablock = np.stack(outs, axis=0)
        pols = np.array(pols).reshape(ntimes, npols, nlays)
        lays = np.array(lays).reshape(ntimes, npols, nlays)
        self.createDimension('time', ntimes)
        if self._vector:
            self.createDimension('points', npoints)
        self.createDimension('layer', nlays)
        self.createDimension('latitude', nlats)
        self.createDimension('longitude', nlons)
        self.createDimension('nv', 2)
        assert((lays[:, [0], :] == lays[:, :, :]).all())
        assert((pols[:, :, [0]] == pols[:, :, :]).all())
        if self._vector:
            poldims = ('points',)
            for key, idx in [('T', 0), ('J', -3), ('I', -2), ('K', -4)]:
                self.createVariable(
                    key, 'i', poldims,
                    units='0-based index',
                    values=datablock[:, idx].astype('i')
                )
        else:
            poldims = ('time', 'layer', 'latitude', 'longitude')
        for pi, pol in enumerate(pols[0, :, 0]):
            if self._vector:
                poldata = datablock[:, -1][datablock[:, 1] == pi]
            else:
                poldata = datablock[:, pi]
            var = self.createVariable(pol.decode(), 'f', poldims,
                                      values=poldata)
            var.units = 'arbitrary'
            var.description = pol.decode()

        time = self.createVariable(
            'time', 'd', ('time',),
            units='seconds since 1970-01-01 00:00:00+0000'
        )
        time_bounds = self.createVariable(
            'time_bounds', 'd', ('time', 'nv'),
            units='seconds since 1970-01-01 00:00:00+0000'
        )
        rdate = strptime('1970-01-01 00:00:00+0000', '%Y-%m-%d %H:%M:%S%z')
        sdt = np.array([(arl2time(YY, MM, DD, HH, mm) - rdate).total_seconds()
                       for YY, MM, DD, HH, mm, FF in starts])
        edt = np.array([(arl2time(YY, MM, DD, HH, mm) - rdate).total_seconds()
                       for YY, MM, DD, HH, mm, FF in stops])
        time_bounds[:, 0] = sdt
        time_bounds[:, 1] = edt
        time[:] = time_bounds.mean(1)


if __name__ == '__main__':
    print(arlconcdump.isMine('cdump24'))
    f = arlconcdump('cdump24')
