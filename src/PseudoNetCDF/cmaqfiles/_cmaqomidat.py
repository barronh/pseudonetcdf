from __future__ import print_function
from ..core._files import PseudoNetCDFFile
import numpy as np
from datetime import datetime, timedelta


class cmaqomidat(PseudoNetCDFFile):
    def __init__(self, path):
        """
        Arguments
        ---------
        path : str
            path to CMAQ-ready OMI file. This is the output of create_omi
            preprocessor. See notes for general form.

        Notes
        -----
        File has the genearl form

        nlat 17
        nlon 17
        yeardate latitude 180Z 157.5W ...
        1990.203       80  426    415 ...
        ...
        1990.203      -80  260    253 ...
        """
        import pandas as pd
        tmpf = open(path, 'r')
        key1, val1 = tmpf.readline().split()
        key2, val2 = tmpf.readline().split()
        if key1.strip() == 'nlat':
            nlat = eval(val1)
            nlon = eval(val2)

        if key1.strip() == 'nlon':
            nlon = eval(val1)
            nlat = eval(val2)

        self._df = df = pd.read_csv(path, delimiter=r'\s+', skiprows=2)

        latitudes = df.latitude.values.reshape(-1, nlat)
        if not (latitudes[0] == latitudes[:]).all():
            raise ValueError('latitudes should all be in the same order')

        # Latitudes are reordered from top-bottom to bottom-top
        data = self._df.drop(
            columns=['yeardate', 'latitude']
        ).values.reshape(-1, nlat, nlon)[:, ::-1]

        ntime = data.shape[0]

        # Latitudes are reordered from top-bottom to bottom-top
        latc = latitudes[0, ::-1]
        dlat = np.abs(np.diff(latc)).mean() / 2
        latb = np.array([latc - dlat, latc + dlat]).T

        lonc = np.linspace(-180, 180, nlon)
        dlon = np.abs(np.diff(lonc)).mean() / 2
        lonb = np.array([lonc - dlon, lonc + dlon]).T

        self.createDimension('time', ntime)
        self.createDimension('latitude', nlat)
        self.createDimension('longitude', nlon)
        self.createDimension('nv', 2)
        fltyears = np.sort(df.yeardate.unique())
        rdate = datetime(1970, 1, 1)
        dts = []
        for fltyear in fltyears:
            year, fyear = divmod(fltyear, 1)
            year = int(year)
            nextyear = datetime(year + 1, 1, 1)
            thisyear = datetime(year, 1, 1)
            ndays = (nextyear - thisyear).total_seconds() / 3600. / 24.
            date = thisyear + timedelta(days=ndays * fyear)
            dt = (date - rdate).total_seconds() / 3600. / 24.
            dts.append(dt)

        dts = np.array(dts)
        self.createVariable(
            'year', 'd', ('time',),
            units='fractional year'
        )[:] = fltyears
        self.createVariable(
            'time', 'd', ('time',),
            units='days since 1970-01-01 00:00:00+0000',
            values=dts
        )
        self.createVariable(
            'latitude', 'd', ('latitude',),
            units='degrees_north', values=latc
        )
        self.createVariable(
            'latitude_bounds', 'd', ('latitude', 'nv'),
            units='degrees_north', values=latb
        )
        self.createVariable(
            'longitude', 'd', ('latitude',),
            units='degrees_east', values=lonc
        )
        self.createVariable(
            'longitude_bounds', 'd', ('latitude', 'nv'),
            units='degrees_east', values=lonb
        )
        v = self.createVariable(
            'O3', 'd', ('time', 'latitude', 'longitude'),
            missing_value=-1, units='DU'
        )
        v[:] = data
        self.setCoords([
            'time', 'latitude', 'longitude',
            'latitude_bounds', 'longitude_bounds'
        ])
