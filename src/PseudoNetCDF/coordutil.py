from PseudoNetCDF import warn
import numpy as np

def gettimes(ifile):
    from datetime import datetime, timedelta
    if 'TFLAG' in ifile.variables.keys():
        dates = ifile.variables['TFLAG'][:][:, 0, 0]
        times = ifile.variables['TFLAG'][:][:, 0, 0]
        yyyy = dates // 1000
        jjj = dates % 1000
        hours = times // 10000
        minutes = times % 10000 // 100
        seconds = times % 100
        days = jjj + (hours + minutes / 60. + seconds / 3600.) / 24.
        out = np.array([datetime(yyyyi, 1, 1) + timedelta(days = day - 1) for day in days])
        return out
    elif 'tau0' in ifile.variables.keys():
        out = datetime(1985, 1, 1, 0) + np.array([timedelta(hours =i) for i in ifile.variables['tau0'][:]])
        return out
    elif 'time' in ifile.variables.keys():
        time = ifile.variables['time']
        if 'since' in time.units:
            unit, base = time.units.strip().split(' since ')
            if 'UTC' in base:
                sdate = datetime.strptime(base, '%Y-%m-%d %H:%M:%S UTC')
            elif 'Z' in base:
                sdate = datetime.strptime(base, '%Y-%m-%d %HZ')
            elif ':' in base:
                sdate = datetime.strptime(base, '%Y-%m-%d %H:%M:%S')
            else:
                sdate = datetime.strptime(base, '%Y-%m-%d')
            out = sdate + np.array([timedelta(**{unit: i}) for i in time[:]])
            return out
        else:
            return time
    else:
        raise ValueError('cannot understand time for file')

def gettimebnds(ifile):
    from datetime import datetime, timedelta
    if 'TFLAG' in ifile.variables.keys():
        dates = ifile.variables['TFLAG'][:][:, 0, 0]
        times = ifile.variables['TFLAG'][:][:, 0, 0]
        yyyy = dates // 1000
        jjj = dates % 1000
        hours = times // 10000
        minutes = times % 10000 // 100
        seconds = times % 100
        days = jjj + (hours + minutes / 60. + seconds / 3600.) / 24.
        out = np.array([datetime(yyyyi, 1, 1) + timedelta(days = day - 1) for day in days])
        
        hours = ifile.TSTEP // 10000
        minutes = ifile.TSTEP % 10000 // 100
        seconds = ifile.TSTEP % 100
        hours = (hours + minutes / 60. + seconds / 3600.)
        return np.array([out, out + timedelta(hours = hours)]).T
    elif 'tau0' in ifile.variables.keys() and 'tau1' in ifile.variables.keys():
        out1 = datetime(1985, 1, 1, 0) + np.array([timedelta(hours =i) for i in ifile.variables['tau0'][:]])
        out2 = datetime(1985, 1, 1, 0) + np.array([timedelta(hours =i) for i in ifile.variables['tau1'][:]])
        return np.array([out1, out2]).T
    elif 'time' in ifile.variables.keys():
        time = ifile.variables['time']
        if 'since' in time.units:
            unit, base = time.units.strip().split(' since ')
            if 'UTC' in base:
                sdate = datetime.strptime(base, '%Y-%m-%d %H:%M:%S UTC')
            elif 'Z' in base:
                sdate = datetime.strptime(base, '%Y-%m-%d %HZ')
            elif ':' in base:
                sdate = datetime.strptime(base, '%Y-%m-%d %H:%M:%S')
            else:
                sdate = datetime.strptime(base, '%Y-%m-%d')
            out = sdate + np.array([timedelta(**{unit: i}) for i in time[:]])
            return np.array([out, out + (out[1] - out[0])]).T
        else:
            return np.array([time, time + (time[1] - time[0])]).T
    else:
        raise ValueError('cannot understand time for file')

def getsigmabnds(ifile):
    if hasattr(ifile, 'VGLVLS'):
        return ifile.VGLVLS[:]
    elif 'etai_pressure' in ifile.variables:
        etai_pressure = ifile.variables['etai_pressure']
        return (etai_pressure - etai_pressure[-1]) / (etai_pressure[0] - etai_pressure[-1])
    elif 'layer_bounds' in ifile.variables:
        lay = ifile.variables['layer_bounds']
        if lay.units.strip() in ('Pa', 'hPa'):
            sigma = (lay[:] -lay[-1]) / (lay[0] - lay[-1])
            return sigma
        else:
            warn("Unknown tranform of layer to sigma; sigma units %s" % lay.units)
            return lay
    else:
        warn("Unknown vertical coordinate")
        if hasattr(ifile, 'NLAYS'):
            nlays = ifile.NLAYS
        elif 'LAY' in ifile.dimensions:
            nlays = len(ifile.dimensions['LAY'])
        elif 'lev' in ifile.dimensions:
            nlays = len(ifile.dimensions['lev'])
        elif 'layer' in ifile.dimensions:
            nlays = len(ifile.dimensions['layer'])
        else:
            nlays = 1
        return np.arange(nlays)
        
def pres_from_sigma(sigma, pref, ptop, avg = False):
    pres = sigma * (pref - ptop) + ptop
    if avg:
        pres = pres[:-1] + np.diff(pres) / 2.
    return pres

def getpresmid(ifile, pref = 101325., ptop = None):
    presb = getpresbnds(ifile, pref = 101325., ptop = None)
    return presb[:-1] + np.diff(presb) / 2

def getsigmamid(ifile):
    sigmab = getsigmabnds(ifile)
    return sigmab[:-1] + np.diff(sigmab) / 2

def getpresbnds(ifile, pref = 101325., ptop = None):
    if 'etai_pressure' in ifile.variables:
        return ifile.variables['etai_pressure'][:]
    elif 'layer_bounds' in ifile.variables:
        return ifile.variables['layer_bounds'][:]
    else:
        sigma = getsigmabnds(ifile)
        if ptop is None:
            if hasattr(ifile, 'VGTOP'):
                ptop = ifile.VGTOP
            else:
                warn("Assuming VGTOP = 10000 Pa")
                ptop = 10000
            
        return pres_from_sigma(sigma, pref = pref, ptop = ptop)

def getlatbnds(ifile):
    if 'latitude_bounds' in ifile.variables:
        latb = ifile.variables['latitude_bounds']
        unit = latb.units.strip()
        if 'nv' in latb.dimensions:
            if latb[:].ndim == 2 and len(ifile.dimensions['nv']) == 2:
                latb = np.append(latb[:][:, 0], latb[:][-1, 1])
            elif latb[:].ndim == 2 and len(ifile.dimensions['nv']) == 4:
                latb = np.append(latb[:][:, 0], latb[:][-1, 1])
            elif latb.ndim == 3:
                latb = latb[:, :, 0]
            
    elif 'latitude' in ifile.variables:
        latb = ifile.variables['latitude']
        unit = latb.units.strip()
        latdiff = np.diff(latb, axis = 0)
        if not (latdiff == latdiff[[0]]).all():
            warn('Latitude bounds are approximate')
        latb = np.concatenate([latb, latb[[-1]]], axis = 0) - .5 * np.concatenate([latdiff[:], latdiff[[-1]], -latdiff[[-1]]], axis = 0)
        if latb.ndim == 2:
            latb = np.append(latb, latb[:, [-1]], axis = 1)
            
    elif 'ROW' in ifile.dimensions:
        unit = 'x (m)'
        latb = np.arange(len(ifile.dimensions['ROW']) + 1) * getattr(ifile, 'YCELL', 1) / 1000.
    else:
        raise KeyError('latitude bounds not found')
    return latb, unit

def getybnds(ifile):
    if 'ROW' in ifile.dimensions:
        unit = 'y (m)'
        latb = np.arange(len(ifile.dimensions['ROW']) + 1) * getattr(ifile, 'YCELL', 1)
    else:
        raise KeyError('latitude bounds not found')
    return latb, unit

def getlonbnds(ifile):
    if 'longitude_bounds' in ifile.variables:
        lonb = ifile.variables['longitude_bounds']
        unit = lonb.units.strip()
        if 'nv' in lonb.dimensions:
            if lonb[:].ndim == 2 and len(ifile.dimensions['nv']) == 2:
                lonb = np.append(lonb[:][:, 0], lonb[:][-1, 1])
            elif lonb[:].ndim == 3:
                lonb = lonb[:][:, :, 0]
    elif 'longitude' in ifile.variables:
        lonb = ifile.variables['longitude']
        unit = lonb.units.strip()
        if lonb.ndim > 1:
            londiff = np.diff(lonb, axis = 1)
            alldiffsame = (londiff == londiff[:, [0]]).all()
        elif lonb.ndim == 1:
            alldiffsame = True
            londiff = np.diff(lonb)
        else:
            raise ValueError("Cannot infer longitude bounds when dimensions >2")
        if not alldiffsame:
            londiff = np.diff(lonb, axis = 1)
            if not (londiff == londiff[:, [0]]).all():
                warn('Longitude bounds are approximate')
            lonb = np.concatenate([lonb, lonb[:, [-1]]], axis = 1) - .5 * np.concatenate([londiff[:, :], londiff[:, [-1]], -londiff[:, [-1]]], axis = 1)
            lonb = np.append(lonb, lonb[[-1], :], axis = 0)
        else:
            londiff = np.diff(lonb, axis = 0)
            lonb = np.concatenate([lonb, lonb[[-1]]], axis = 0) - .5 * np.concatenate([londiff[:], londiff[[-1]], -londiff[[-1]]], axis = 0)

    else:
        raise KeyError('longitude bounds not found')
    return lonb, unit

def getxbnds(ifile):
    if 'COL' in ifile.dimensions:
        unit = 'x (m)'
        lonb = np.arange(len(ifile.dimensions['COL']) + 1) * getattr(ifile, 'XCELL', 1)
    else:
        raise KeyError('x bounds not found')
    return lonb, unit

def getmap(ifile, resolution = 'i'):
    from mpl_toolkits.basemap import Basemap
    from conventions.ioapi import get_ioapi_sphere
    if getattr(ifile, 'GDTYP', 0) in (2, 7) and all([hasattr(ifile, k) for k in 'P_GAM P_ALP P_BET XORIG YORIG XCELL YCELL'.split()]):
        try:
            NROWS = len(ifile.dimensions['ROW'])
            NCOLS = len(ifile.dimensions['COL'])
        except KeyError:
            NROWS = ifile.NROWS
            NCOLS = ifile.NCOLS
            
        llcrnrx = ifile.XORIG
        urcrnrx = ifile.XORIG + NCOLS * ifile.XCELL

        llcrnry = ifile.YORIG
        urcrnry = ifile.YORIG + NROWS * ifile.YCELL
        semi_major_axis, semi_minor_axis = get_ioapi_sphere()
        if ifile.GDTYP == 2:
            from mpl_toolkits.basemap.pyproj import Proj
            p = Proj(proj='lcc',rsphere = (semi_major_axis, semi_major_axis), lon_0 = ifile.P_GAM, lat_1 = ifile.P_ALP, lat_2 = ifile.P_BET, lat_0 = ifile.YCENT)
            llcrnrlon, llcrnrlat = p(llcrnrx, llcrnry, inverse = True)
            urcrnrlon, urcrnrlat = p(urcrnrx, urcrnry, inverse = True)
            m = Basemap(projection = 'lcc', rsphere = (semi_major_axis, semi_major_axis), lon_0=ifile.P_GAM, lat_1 = ifile.P_ALP, lat_2 = ifile.P_BET, lat_0 = ifile.YCENT, llcrnrlon = llcrnrlon, llcrnrlat = llcrnrlat, urcrnrlat = urcrnrlat, urcrnrlon = urcrnrlon, resolution = resolution, suppress_ticks = False)
        elif ifile.GDTYP == 7:
            from mpl_toolkits.basemap.pyproj import Proj
            mapstr = '+proj=merc +a=%s +b=%s +lat_ts=0 +lon_0=%s' % (semi_major_axis, semi_major_axis, ifile.XCENT)
            p = Proj(mapstr)
            #p = Proj(proj='merc',rsphere = (semi_major_axis, semi_major_axis), lat_ts = ifile.P_ALP, lat_0 = ifile.YCENT, lon_0 = ifile.XCENT)
            llcrnrlon, llcrnrlat = p(llcrnrx, llcrnry, inverse = True)
            urcrnrlon, urcrnrlat = p(urcrnrx, urcrnry, inverse = True)
            m = Basemap(projection = 'merc', rsphere = (semi_major_axis, semi_major_axis), lon_0=ifile.XCENT, lat_ts = 0, llcrnrlon = llcrnrlon, llcrnrlat = llcrnrlat, urcrnrlat = urcrnrlat, urcrnrlon = urcrnrlon, resolution = resolution, suppress_ticks = False)
        print 'Found IO/API Mapping parameters'
    else:
        kwds = dict(suppress_ticks = False)
        try:
            lat, latunit = getlatbnds(ifile)
            lon, lonunit = getlonbnds(ifile)
            kwds['resolution'] = 'i'
            kwds['llcrnrlat'] = lat[:].min()
            kwds['urcrnrlat'] = lat[:].max()
            kwds['llcrnrlon'] = lon[:].min()
            kwds['urcrnrlon'] = lon[:].max()
            kwds['resolution'] = resolution
        except:
            pass
        m = Basemap(**kwds)
    return m

