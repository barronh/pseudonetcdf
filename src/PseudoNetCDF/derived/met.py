from PseudoNetCDF.sci_var import PseudoNetCDFMaskedVariable as PseudoNetCDFVariable
import numpy as np

def add_derived_met(ifile, func, *args, **kwds):
    results = func(*args, **kwds)
    if not isinstance(results, tuple):
        results = (results,)
    for result in results:
        ifile.variables[result.short_name] = result

def wmr_ptd( p, td, gkg = False):
    """
    Inputs:
       p = surface pressure in mb; 
       Td = dew point in deg C; 
    Calculates
       e = vapor pressure in mb; 

    Returns
       q = specific humidity in kg/kg. 
    """
    dimensions = getattr(p, 'dimensions', getattr(td, 'dimensions', 'unknown'))
    e = 6.112*np.ma.exp((17.67*td)/(td + 243.5)); 
    q = (0.622 * e)/(p - (0.378 * e));
    if gkg:
        q *= 1000.
        units = 'g/kg'
    else:
        units = 'kg/kg'
    Q = PseudoNetCDFVariable(None, 'Q', q.dtype.char, dimensions, units = units, long_name = 'specific humidity in kg/kg', short_name = 'Q', values = q)
    return Q

def relhum_ttd(t, td, percent = False, add = True):
    """
    t - temperature in K
    td - dewpoint temperature in K


    ;
    ; Calculate relative humidity given temperature (K)
    ; and dew point temperature (K)
    ;
    ; reference: John Dutton, Ceaseless Wind, 1976
    """      
    dimensions = getattr(t, 'dimensions', getattr(td, 'dimensions', ('unknown',)))
    gc  = 461.5                        # [j/{kg-k}]   gas constant water vapor
    gc  = gc/(1000.*4.186)             # [cal/{g-k}]  change units
                                       # lhv=latent heat vap
    lhv = ( 597.3-0.57*(t-273.) )      # dutton top p273 [empirical]

    values = np.ma.exp( (lhv/gc)*(1.0/t - 1.0/td))
    if percent:
        values *= 100.
    
    rh           = PseudoNetCDFVariable(None, 'RH', 'f', dimensions, values = values )
    rh.long_name = "relative humidity"
    if percent:
        rh.units = "%"
    else:
        rh.units     = "fraction"
    rh.short_name = 'RH'
    return rh

def wmrq_ptd(p, td, dry = False, kgkg = False):
    """
    p - pressure in Pa
    td - temperature in degrees celcius
    dry - return results in mixing ratio or specific humidity in dry air
    kgkg - return results in kg/kg instead of g/kg
    """
    dimensions = getattr(p, 'dimensions', getattr(td, 'dimensions', ('unknown',)))
    p = np.asarray(p)
    td = np.asarray(td)
    # ncl:   q = mixhum_ptd (p,td,option)   [ q=g/kg ]
    #                        p  - PA
    #                        td - K
    # local
    # INTEGER  N  
    # DOUBLE PRECISION T0, PA2MB
    # DATA     T0    /273.15d0/
    # DATA     PA2MB /0.01d0  /
    T0 = 273.15
    PA2MB = .01

    # mixing ratio (kg/kg)
    # the function wants hPA (mb) and degrees centigrade
    WMR = wmr_skewt_pt(p*PA2MB, td - T0) * 0.001

    # if iswit=2 calculate specific humidity (kg/kg)
      
    if dry:
        name = 'W'
        long_name = 'mass of water per mass of dry air'
    else:
        name = 'Q'
        long_name = 'mass of water per mass of air (dry+wet)'
        WMR = WMR/(WMR + 1.0)

    # if iswit < 0 then return g/kg
      
    if kgkg:
        WMR *= 1000.
        units = 'g/kg'
    else:
        units = 'kg/kg'
    wmr           = PseudoNetCDFVariable(None, name, 'f', dimensions, values = WMR, units = units, long_name = long_name )
    wmr.short_name = name
    return wmr

def wmr_skewt_pt(p,t):
    """
    p - pressure in mb
    t - temperature in degrees C
    this function approximates the mixing ratio wmr (grams of water
    vapor per kilogram of dry air) given the pressure p (mb) and the
    temperature t (celsius). the formula used is given on p. 302 of the
    smithsonian meteorological tables by roland list (6th edition).
    """
    #eps = ratio of the mean molecular weight of water (18.016 g/mole)
    #      to that of dry air (28.966 g/mole)
    EPS = 0.62197

    #the next two lines contain a formula by herman wobus for the
    #correction factor wfw for the departure of the mixture of air
    #and water vapor from the ideal gas law. the formula fits values
    #in table 89, p. 340 of the smithsonian meteorological tables,
    #but only for temperatures and pressures normally encountered in
    #in the atmosphere.

    X = 0.02* (t-12.5+7500./p)
    WFW = 1. + 4.5e-06*p + 1.4e-03*X*X
    FWESW = WFW*esw_skewt_t(t)
    R = EPS*FWESW/ (p-FWESW)

    #convert r from a dimensionless ratio to grams/kilogram.

    return 1000.*R

def esw_skewt_t(t):
    """
    Arguments:
        t - temperature in degrees C
    Returns:
        esw - Saturation vapor pressure in millibars
    Adapted from NCL
    this function returns the saturation vapor pressure esw (millibars)
    over liquid water given the temperature t (celsius). the polynomial
    approximation below is due to herman wobus, a mathematician who
    worked at the navy weather research facility, norfolk, virginia,
    but who is now retired. the coefficients of the polynomial were
    chosen to fit the values in table 94 on pp. 351-353 of the smith-
    sonian meteorological tables by roland list (6th edition). the
    approximation is valid for -50 < t < 100c.
 
    es0 = saturation vapor ressure over liquid water at 0c
    """
    ES0 = 6.1078e0

    POL = 0.99999683e0 + t* (-0.90826951e-02+
          t* (0.78736169e-04+t* (-0.61117958e-06+t* (0.43884187e-08+
          t* (-0.29883885e-10+t* (0.21874425e-12+t* (-0.17892321e-14+
          t* (0.11112018e-16+t* (-0.30994571e-19)))))))))
    return ES0/POL**8

def uv_wswd(ws, wd, isradians = False):
    """
    Arguments:
        ws - wind speed array
        wd - wind direction array (in degrees unless isradians is set to true)
        isradians - False if WD is in degrees
    Returns:
        U, V - u-component and v-component winds
    
    https://www.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html
    """
    dimensions = getattr(u, 'dimensions', getattr(v, 'dimensions', ('unknown',)))
    units = getattr(ws, 'units', 'unknown')
    if isradians:
        direction = wd
    else:
        direction = np.radians(wd)
    wind_speed = ws
    U = -np.sin(direction) * wind_speed; 
    V = -np.cos(direction) * wind_speed;
    u = PseudoNetCDFVariable(None, 'U', U.dtype.char, dimensions, values = U, units = units, short_name = 'U')
    v = PseudoNetCDFVariable(None, 'V', V.dtype.char, dimensions, values = V, units = units, short_name = 'V')
    return u, v
    
def wswd_uv(u, v, return_radians = False):
    """
    Arguments:
        u - array of u-component winds
        v - array of v-domponent winds
        return_radians - convert direction to radians
    Returns:
        ws, wd - wind speed and direction (from direction)
        
    https://www.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html
    """
    dimensions = getattr(u, 'dimensions', getattr(v, 'dimensions', ('unknown',)))
    units = getattr(u, 'units', getattr(v, 'units', 'unknown'))
    wind_speed = np.sqrt(u**2 + v**2); 
    uvarctan = np.arctan(u/v)
    wind_direction = np.where(v<0, uvarctan*180./np.pi, uvarctan * 180./np.pi + 180.)
    if return_radians:
        wind_direction = np.radians(wind_direction)
    ws = PseudoNetCDFVariable(None, 'WS', wind_speed.dtype.char, dimensions, values = wind_speed, units = units, short_name = 'WS')
    wd = PseudoNetCDFVariable(None, 'WD', wind_speed.dtype.char, dimensions, values = wind_direction, units = 'radians' if return_radians else 'degrees', short_name = 'WD')
    return ws, wd