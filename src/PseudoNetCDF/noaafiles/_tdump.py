from PseudoNetCDF import PseudoNetCDFFile
import numpy as np

_units = dict(trajid = '---',
              metgridid = '---',
              year = 'year',
              month = 'month of year',
              day = 'day of month',
              hour = 'hour of day',
              minute = 'minute of hour',
              forecast_hour = 'hour in forecast',
              age = 'hours',
              pressure = 'hPa',
              theta = 'K',
              air_temp = 'K',
              rainfall = 'mm/h',
              mixdepth = 'm',
              relhumid = '%',
              terr_msl = 'm',
              sun_flux = 'W/m**2',)

def _timefromnoaa(year, month, day, hour, minute):
    from datetime import datetime, timezone
    datei = (year.astype('l') * 100000000 +
             month.astype('l') * 1000000 +
             day.astype('l') * 10000 +
             hour.astype('l') * 100 +
             minute.astype('l'))
    datestrs = np.char.decode(np.char.add(datei.astype('S16'), b'+0000'))
    dates = np.array([datetime.strptime(d, '%y%m%d%H%M%z') for d in datestrs])
    return dates

def _year(year):
    if np.floor(np.log10(year)) == 1.:
       if year < 70:
           return year + 1900
       else:
           return year + 2000
    else:
       return year

_vyear = np.vectorize(_year)

class arltrajdump(PseudoNetCDFFile):
    @classmethod
    def isMine(cls, path):
        try:
            f = arltrajdump(path)
            return True
        except:
            return False
        
    def __init__(self, path):
        self._path = path
        f = self._file = open(path)
        """
        Record #1

        I6 - Number of meteorological grids used in calculation
        """
        nmgrid = self.NMETGRIDS = int(f.readline()[:6].strip())
        """
        Records Loop #2 through the number of grids

        A8 - Meteorological Model identification
        5I6 - Data file starting Year, Month, Day, Hour, Forecast Hour
        """
        metgridlines = [f.readline().strip().split() for i in range(nmgrid)]
        self.METGRIDS = ','.join([l[0] for l in metgridlines])
        self.createDimension('metgrid', nmgrid)
        v = self.createVariable('met_year', 'i', ('metgrid',))
        v.units = 'year'
        v.long_name = 'year'
        v[:] = [int(l[1]) for l in metgridlines]
        v = self.createVariable('met_month', 'i', ('metgrid',))
        v.units = 'month'
        v.long_name = 'month of the year'
        v[:] = [int(l[2]) for l in metgridlines]
        v = self.createVariable('met_day', 'i', ('metgrid',))
        v.units = 'day'
        v.long_name = 'day of the month'
        v[:] = [int(l[3]) for l in metgridlines]
        v = self.createVariable('met_hour', 'i', ('metgrid',))
        v.units = 'hour'
        v.long_name = 'hour of the day (GMT)'
        v[:] = [int(l[4]) for l in metgridlines]
        v = self.createVariable('met_forecast_hour', 'i', ('metgrid',))
        v.units = 'forecast_hour'
        v.long_name = 'hour of the forecast'
        v[:] = [int(l[5]) for l in metgridlines]
        
        """
        Record #3

        I6 - number of different trajectories in file
        1X,A8 - direction of trajectory calculation (FORWARD, BACKWARD)
        1X,A8 - vertical motion calculation method (OMEGA, THETA, ...)
        """
        ntrajstr, self.TRAJECTORY_CALCULATION, self.VERTMOTION = f.readline().strip().split()
        ntrajs = self.NTRAJECTORIES = int(ntrajstr)
        """
        Record Loop #4 through the number of different trajectories in file

        4I6 - starting year, month, day, hour
        2F9.3 - starting latitude, longitude
        F8.1 - starting level above ground (meters)
        """
        trajmeta = np.array([f.readline().strip().split() for i in range(ntrajs)], dtype = 'f')
        self.createDimension('trajectory', ntrajs)
        v = self.createVariable('traj_year', 'i', ('trajectory',))
        v.units = 'year'
        v.long_name = 'year'
        v[:] = trajmeta[:, 0].astype('i')
        v = self.createVariable('traj_month', 'i', ('trajectory',))
        v.units = 'month'
        v.long_name = 'month of the year'
        v[:] = trajmeta[:, 1].astype('i')
        v = self.createVariable('trajectory_day', 'i', ('trajectory',))
        v.units = 'day'
        v.long_name = 'day of the month'
        v[:] = trajmeta[:, 2].astype('i')
        v = self.createVariable('trajectory_hour', 'i', ('trajectory',))
        v.units = 'hour'
        v.long_name = 'hour of the day (GMT)'
        v[:] = trajmeta[:, 3].astype('i')
        v = self.createVariable('trajectory_init_latitude', 'i', ('trajectory',))
        v.units = 'degrees_north'
        v.long_name = 'initial latitude'
        v[:] = trajmeta[:, 4]
        v = self.createVariable('trajectory_init_longitude', 'i', ('trajectory',))
        v.units = 'degrees_east'
        v.long_name = 'initial longitude'
        v[:] = trajmeta[:, 5]
        v = self.createVariable('trajectory_init_height', 'i', ('trajectory',))
        v.units = 'meters agl'
        v.long_name = 'initial altitude'
        v[:] = trajmeta[:, 6]
        # Starting time
        self._starttimes = _timefromnoaa(trajmeta[:, 0], trajmeta[:, 1], trajmeta[:, 2], trajmeta[:, 3], trajmeta[:, 3]*0)
        """
        Record #5

        I6 - number (n) of diagnostic output variables
        n(1X,A8) - label identification of each variable (PRESSURE, THETA, ...)
        """
        diagnostics = f.readline().strip().lower().split()
        ndiag = int(diagnostics[0])
        assert(ndiag+1 == len(diagnostics))
        """
        Record Loop #6 through the number of hours in the simulation

        I6 - trajectory number
        I6 - meteorological grid number or antecedent trajectory number
        5I6 - year month day hour minute of the point
        I6 - forecast hour at point
        F8.1 - age of the trajectory in hours
        2F9.3 - position latitude and longitude
        1X,F8.1 - position height in meters above ground
        n(1X,F8.1) - n diagnostic output variables (1st to be output is always pressure)
        """
        try:
            import pandas as pd
        except:
            raise ImportError('arltrajdump requires pandas; install pandas (e.g., pip install pandas)')
        data = pd.read_csv(f, delimiter = '\s+', names = 'trajid metgridid year month day hour minute forecast_hour age latitude longitude altitude'.split() + diagnostics[1:])#, parse_dates = ['YEAR MONTH DAY HOUR MINUTE'.split()])
        mytimes = _timefromnoaa(data['year'], data['month'], data['day'],
                                data['hour'], data['minute'])
        unique_times = np.sort(np.unique(mytimes))
        ntimes = len(unique_times)
        self.createDimension('time', ntimes)
        utraj = data['trajid'].unique()
        mytraj = data['trajid'].values
        myage = data['age'].values
        trajidx = (utraj[:, None] == mytraj[None,:]).argmax(0)
        timeidx = (unique_times[:, None] == mytimes[None,:]).argmax(0)
        
        tmpv = np.ma.masked_all((ntimes, ntrajs), dtype = 'f')
        for k in data.columns:
            v = self.createVariable(k, 'f', ('time', 'trajectory'), fill_value = -999.)
            v.long_name = k
            v.units = _units.get(k, 'unknown')
            v[:] = tmpv
            v[timeidx, trajidx] = data[k].values
        #self._data = data
        
    def getTimes(self):
        year = (self.variables['year']).max(1).astype('l')
        month = self.variables['month'].max(1).ravel().astype('l')
        day = self.variables['day'].max(1).astype('l')
        hour = self.variables['hour'].max(1).astype('l')
        minute = self.variables['minute'].max(1).astype('l')
        return _timefromnoaa(year, month, day, hour, minute)

if __name__ == '__main__':
    f = arltrajdump('tdump_008')
