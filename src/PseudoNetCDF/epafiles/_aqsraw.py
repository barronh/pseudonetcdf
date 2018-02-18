from __future__ import print_function
from PseudoNetCDF.sci_var import PseudoNetCDFFile
import numpy as np
from datetime import datetime, timedelta

def getbdate(x):
    if x is None:
        return x
    elif isinstance(x, (datetime, np.datetime64)):
        return x
    else:
        return datetime.strptime(x + ' 00:00', '%Y-%m-%d %H:%M')

def getedate(x):
    if x is None:
        return x
    elif isinstance(x, (datetime, np.datetime64)):
        return x
    else:
        return datetime.strptime(x + ' 23:59', '%Y-%m-%d %H:%M')


class aqsraw(PseudoNetCDFFile):
    def __init__(self, yearpath, timeresolution = 'hourly', bdate = None, edate = None, rdate = datetime(1900, 1, 1), wktpolygon = None, sampleval = None, verbose = 0):
        """
        yearpath - path to csv file from AQS
        timeresolution - choices = ['daily', 'hourly'] default = 'hourly', Defaults to hourly'
        bdate - Limit output so that the date starts at YYYY-MM-DD
        edate - Limit output so that the date end YYYY-MM-DD (inclusive)
        rdate - datetime(1900, 1, 1), dest = 'rdate', type = getrdate, help = 'Reference date YYYYMMDD HH:MM:SS')
        wktpolygon - WKT Polygon (default: None) equivalent to "POLYGON ((-180 -90, 180 -90, 180 90, -180 90, -180 -90))"
        sampleval - Defaults to "Sample Measurement" for hourly and "Arithmetic Mean" for daily
        verbose - level of verbosity
        """
        try:
            import pandas
        except:
            raise ImportError('aqsraw requires pandas; install pandas (e.g., pip install pandas)')
        if not wktpolygon is None:
            from shapely.wkt import loads
            from shapely.prepared import prep
            bounds = loads(wktpolygon)
            pbounds = prep(bounds)
        else:
            from shapely.geometry import Point
        bdate = getbdate(bdate)
        edate = getedate(edate)
        nseconds = {'hourly': 3600, 'daily': 3600*24}[timeresolution]
        tunit = {'hourly': 'hours', 'daily': 'days'}[timeresolution]

        if sampleval is None:
            if timeresolution == 'hourly':
                sampleval = "Sample Measurement"
            elif timeresolution == 'daily':
                sampleval = "Arithmetic Mean"
            else:
                raise KeyError(sampleval + ' not appropriate sampleval value')
        hourlys = []
        for yearpath in [yearpath]:
            if verbose > 0: print('Reading', yearpath)
            if timeresolution == 'hourly':
                parse_dates = [[11, 12]]
                date_key = 'Date GMT_Time GMT'
            else:
                parse_dates = [11]
                date_key = 'Date Local'
            data = pandas.read_csv(yearpath, index_col = False, converters = {'State Code': str, 'County Code': str, 'Site Num': str}, parse_dates = parse_dates)
            intimes = np.array([True]).repeat(len(data), 0)
            #intimes = intimes & (data['Parameter Code'] == param)
            
            if not bdate is None:
                intimes = intimes & (data[date_key] >= bdate)
            if not edate is None:
                intimes = intimes & (data[date_key] < edate)
            
            if wktpolygon is None:
                inspace_and_time = intimes
            else:
                inspace_and_time = np.array([False]).repeat(len(data), 0)
                for idx, (intime, plon, plat) in enumerate(zip(intimes, data['Longitude'].values, data['Latitude'].values)):
                    if intime:
                        inspace_and_time[idx] = pbounds.contains(Point(plon, plat))
            mask = inspace_and_time
       

            hourly = data[mask].groupby(['Parameter Code', 'Parameter Name', 'Units of Measure', date_key, 'State Code', 'County Code', 'Site Num'], as_index = False).aggregate(np.mean).sort_values(by = ['Parameter Code', 'Parameter Name', 'Units of Measure', date_key, 'State Code', 'County Code', 'Site Num'])
    
            hourlys.append(hourly)

        if verbose > 0: print('Concatenating files')
        if len(hourlys) > 1:
            hourly = pandas.concat(hourlys)
        else:
            hourly = hourlys[0]


        sites = hourly.groupby(['State Code', 'County Code', 'Site Num'], as_index = False).aggregate(np.mean)
        nsites = len(sites)

        rawdates = hourly[date_key]
        if bdate is None:
            bdate = rawdates.min()
        if edate is None:
            edate = rawdates.max()
        ntimes = int((edate - bdate).total_seconds() // nseconds) + 1
        alltimes = [timedelta(**{tunit: i}) + bdate for i in range(ntimes)]

        if verbose > 0: print('Creating output file')
        tdim = self.createDimension('time', ntimes)
        tdim.setunlimited(True)
        self.createDimension('LAY', 1)
        sitelist = [row['State Code'] + row['County Code'] + row['Site Num'] for idx, row in sites.iterrows()]
        self.SITENAMES = ';'.join(sitelist)

        self.createDimension('points', len(sitelist))
        self.lonlatcoords = '/'.join(['%f,%f' % (row['Longitude'], row['Latitude']) for idx, row in sites.iterrows()])
        last_time = start_time = hourly[date_key].values.min()
        end_time = hourly[date_key].values.max()
        temp = {}
        lat = self.createVariable('latitude', 'f', ('points',))
        lat.units = 'degrees_north'
        lat.standard_name = 'latitude'
        lat[:] = sites['Latitude'].values
        lon = self.createVariable('longitude', 'f', ('points',))
        lon.units = 'degrees_east'
        lon.standard_name = 'longitude'
        lon[:] = sites['Longitude'].values

        if verbose > 0: print('Processing data rows')

        for idx, row in hourly.iterrows():
            this_time = row[date_key]
            val = row[sampleval]
            unit = row['Units of Measure'].strip()
            aqsid = row['State Code'] + row['County Code'] + row['Site Num']
            sidx = sitelist.index(aqsid)
            var_name = row['Parameter Name']
            if not var_name in self.variables.keys():
                var = self.createVariable(var_name, 'f', ('time', 'LAY', 'points'), fill_value = -999)
                var.units = unit
                var.standard_name = var_name
                temp[var_name] = np.ma.masked_all((ntimes, 1, nsites), dtype = 'f')
                temp[var_name].set_fill_value(-999)
            tmpvar = temp[var_name]
            var = self.variables[var_name]
            assert(var.units == unit)
            if last_time != this_time:
                last_time = this_time
                tidx = alltimes.index(this_time.to_pydatetime())
                if verbose > 0: print('\b\b\b\b\b\b %3d%%' % int(tidx / float(ntimes) * 100), end = '', flush = True)
            if tidx > ntimes:
                raise ValueError('Times (%d) exceed expected (%d)' % (tidx, ntimes))
    
            tmpvar[tidx, 0, sidx] = val
        if verbose > 0: print('')
        if verbose > 0: print('Writing to file')
        for varkey, tempvals in temp.items():
            self.variables[varkey][:] = tempvals

        outtimes = np.array([(t - rdate).total_seconds() / nseconds for t in alltimes])
        time = self.createVariable('time', outtimes.dtype.char, ('time',))
        time.units = rdate.strftime(tunit + ' since %F')
        time.standard_name = 'time'
        time[:] = outtimes

if __name__ == '__main__':
    import sys
    f = aqsraw(sys.argv[1], bdate = datetime(2015,1,1), edate = datetime(2015,2,1))
