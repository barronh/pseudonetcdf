from warnings import warn
import numpy as np

def getsigmabnds(ifile):
    if hasattr(ifile, 'VGLVLS'):
        return ifile.VGLVLS[:]
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
            nlays = len(ifile.dimensions['laer'])
        else:
            nlays = 1
        return np.arange(nlays)
        
def getpresmid(ifile):
    presb = getpresbnds(ifile)
    return presb[:-1] + np.diff(presb) / 2

def getsigmamid(ifile):
    sigmab = getsigmabnds(ifile)
    return sigmab[:-1] + np.diff(sigmab) / 2

def getpresbnds(ifile):
    if 'layer_bounds' in ifile.variables:
        return ifile.variables['layer_bounds'][:]
    else:
        sigma = getsigmabnds(ifile)
        if hasattr(ifile, 'VGTOP'):
            ptop = ifile.VGTOP
        else:
            warn("Assuming VGTOP = 100 hPa")
            ptop = 10000
        return pres_from_sigma(sigma, pref = 101325., ptop = ptop)

def getlatbnds(ifile):
    if 'latitude_bounds' in ifile.variables:
        latb = ifile.variables['latitude_bounds']
        unit = latb.units.strip()
        if 'nv' in latb.dimensions and latb[:].ndim == 2:
            if len(ifile.dimensions['nv']) == 2:
                latb = np.append(latb[:][:, 0], latb[:][-1, 1])
            elif len(ifile.dimensions['nv']) == 4:
                latb = np.append(latb[:][:, 0], latb[:][-1, 1])
            elif latb.ndim == 3:
                latb = latb[:, :, 0]
            
    elif 'latitude' in ifile.variables:
        latb = ifile.variables['latitude']
        unit = latb.units.strip()
        latdiff = np.diff(latb, axis = 0)
        if (latdiff == latdiff[0]).all():
            latb = latb - .5 * latdiff[0]
            latb = np.append(latb, latb[-1] + latdiff[0])
    elif 'ROW' in ifile.dimensions:
        unit = 'x (LCC m)'
        latb = np.arange(len(ifile.dimensions['ROW']) + 1) * getattr(ifile, 'YCELL', 1) / 1000.
    else:
        raise KeyError('latitude bounds not found')
    return latb, unit

def getybnds(ifile):
    if 'ROW' in ifile.dimensions:
        unit = 'y (LCC m)'
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
        londiff = np.diff(lonb, axis = 0)
        if (londiff == londiff[0]).all():
            lonb = lonb - .5 * londiff[0]
            lonb = np.append(lonb, lonb[-1] + londiff[0])
    else:
        raise KeyError('longitude bounds not found')
    return lonb, unit

def getxbnds(ifile):
    if 'COL' in ifile.dimensions:
        unit = 'x (LCC m)'
        lonb = np.arange(len(ifile.dimensions['COL']) + 1) * getattr(ifile, 'XCELL', 1)
    else:
        raise KeyError('x bounds not found')
    return lonb, unit

def getmap(ifile):
    from mpl_toolkits.basemap import Basemap
    from conventions.ioapi import get_ioapi_sphere
    if getattr(ifile, 'GDTYP', 0) != 1 and all([hasattr(ifile, k) for k in 'P_GAM P_ALP P_BET XORIG YORIG XCELL YCELL'.split()]):
        llcrnrx = ifile.XORIG
        urcrnrx = ifile.XORIG + len(ifile.dimensions['COL']) * ifile.XCELL

        llcrnry = ifile.YORIG
        urcrnry = ifile.YORIG + len(ifile.dimensions['ROW']) * ifile.YCELL
        semi_major_axis, semi_minor_axis = get_ioapi_sphere()
        if ifile.GDTYP == 2:
            m = Basemap(projection = 'lcc', a = semi_major_axis, b = semi_major_axis, lon_0=ifile.P_GAM, lat_1 = ifile.P_ALP, lat_2 = ifile.P_BET, lat_0 = ifile.YCENT, llcrnrx = llcrnrx, llcrnry = llcrnry, urcrnry = urcrnry, urcrnrx = urcrnrx, resolution = 'i', suppress_ticks = False)
        
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
        except:
            pass
        m = Basemap(**kwds)
    return m

