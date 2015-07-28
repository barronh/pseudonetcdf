__all__ = ['NO', 'NP', 'NOP', 'MO', 'MP', 'MdnO', 'MdnP', 'STDO', 'STDP', 'RM', 'RMdn', 'MB', 'MdnB', 'WDMB', 'WDMdnB', 'FB', 'MNB', 'MdnNB', 'NMB', 'NMdnB', 'USUTPB', 'PSUTMNPB', 'PSUTMdnNPB', 'PSUTNMPB', 'PSUTNMdnPB', 'ME', 'MdnE', 'WDME', 'WDMdnE', 'FE', 'MNE', 'MdnNE', 'NME', 'NMdnE', 'USUTPE',  'PSUTMNPE', 'PSUTMdnNPE', 'PSUTNMPE', 'PSUTNMdnPE', 'R2', 'RMSE', 'RMSEs', 'RMSEu', 'E1', 'IOA', 'd1', 'AC', 'WDIOA', 'WDRMSE', 'WDAC']
from pncload import createconsole
import numpy as np

def STDO(obs, mod, axis = None):
    """ Standard deviation of Observations """
    return np.ma.std(obs, axis = axis)

def STDP(obs, mod, axis = None):
    """ Standard deviation of Predictions """
    return np.ma.std(mod, axis = axis)

def MNB(obs, mod, axis = None):
    """ Mean Normalized Bias (%)"""
    return  np.ma.masked_invalid((mod-obs)/obs).mean(axis = axis)*100.

def MNE(obs, mod, axis = None):
    """ Mean Normalized Gross Error (%)"""
    return np.ma.masked_invalid(np.ma.abs(mod-obs)/obs).mean(axis = axis)*100.

def MdnNB(obs, mod, axis = None):
    """ Median Normalized Bias (%)"""
    return np.ma.median(np.ma.masked_invalid((mod-obs)/obs), axis = axis)*100.

def MdnNE(obs, mod, axis = None):
    """ Median Normalized Gross Error (%)"""
    return np.ma.median(np.ma.masked_invalid(np.ma.abs(mod-obs)/obs), axis = axis)*100.

def NMdnE(obs, mod, axis = None):
    """ Normalized Median Gross Error (%)"""
    return np.ma.masked_invalid(np.ma.abs(mod-obs).mean(axis = axis)/obs.mean(axis = axis))*100.

def NO(obs, mod, axis = None):
    """ N Observations (#)"""
    return (np.ma.getmaskarray(obs) == False).sum(axis = axis)

def NOP(obs, mod, axis = None):
    """ N Observations/Prediction Pairs (#)"""
    obsc, modc = matchmasks(obs, mod)
    return (np.ma.getmaskarray(obsc) == False).sum(axis = axis)

def NP(obs, mod, axis = None):
    """ N Predictions (#)"""
    return (np.ma.getmaskarray(mod) == False).sum(axis = axis)


def MO(obs, mod, axis = None):
    """ Mean Observations (obs unit)"""
    return obs.mean(axis = axis)

def MP(obs, mod, axis = None):
    """ Mean Predictions (model unit)"""
    return mod.mean(axis = axis)

def MdnO(obs, mod, axis = None):
    """ Median Observations (obs unit)"""
    return np.ma.median(obs, axis = axis)

def MdnP(obs, mod, axis = None):
    """ Median Predictions (model unit)"""
    return np.ma.median(mod, axis = axis)

def RM(obs, mod, axis = None):
    """ Mean Ratio Observations/Predictions (none)"""
    return np.ma.masked_invalid(obs/mod).mean(axis = axis)

def RMdn(obs, mod, axis = None):
    """ Median Ratio Observations/Predictions (none)"""
    return np.ma.median(np.ma.masked_invalid(obs/mod), axis = axis)

def MB(obs, mod, axis = None):
    """ Mean Bias"""
    return (mod-obs).mean(axis = axis)

def MdnB(obs, mod, axis = None):
    """ Median Bias"""
    return np.ma.median(mod-obs, axis = axis)

def WDMB(obs, mod, axis = None):
    """ Wind Direction Mean Bias"""
    return circlebias(mod-obs).mean(axis = axis)

def WDMdnB(obs, mod, axis = None):
    """ Wind Direction Median Bias"""
    return np.ma.median(circlebias(mod-obs), axis = axis)

def NMB(obs, mod, axis = None):
    """ Normalized Mean Bias (%)"""
    return (mod-obs).sum(axis = axis)/obs.sum(axis = axis)*100.

def NMdnB(obs, mod, axis = None):
    """ Normalized Median Bias (%)"""
    return np.ma.median(mod-obs, axis = axis)/np.ma.median(obs, axis = axis)*100.

def FB(obs, mod, axis = None):
    """ Fractional Bias (%)"""
    return ((np.ma.masked_invalid((mod-obs)/(mod+obs))).mean(axis = axis)*2.)*100.

def ME(obs, mod, axis = None):
    """ Mean Gross Error (model and obs unit)"""
    return np.ma.abs(mod-obs).mean(axis = axis)

def MdnE(obs, mod, axis = None):
    """ Median Gross Error (model and obs unit)"""
    return np.ma.median(np.ma.abs(mod-obs), axis = axis)

def WDME(obs, mod, axis = None):
    """ Wind Direction Mean Gross Error (model and obs unit)"""
    return np.ma.abs(circlebias(mod-obs)).mean(axis = axis)

def WDMdnE(obs, mod, axis = None):
    """ Wind Direction Median Gross Error (model and obs unit)"""
    cb = circlebias(mod-obs)
    return np.ma.median(np.ma.abs(cb), axis = axis)

def NME(obs, mod, axis = None):
    """ Normalized Mean Error (%)"""
    out = (np.ma.abs(mod-obs).sum(axis = axis)/obs.sum(axis = axis))*100
    return out

def NMdnE(obs, mod, axis = None):
    """ Normalized Median Error (%)"""
    out = np.ma.median(np.ma.abs(mod-obs), axis = axis)/np.ma.median(obs, axis = axis)*100
    return out

def FE(obs, mod, axis = None):
    """ Fractional Error (%)"""
    return (np.ma.abs(mod-obs)/(mod+obs)).mean(axis = axis)*2.*100.

def USUTPB(obs, mod, axis = None):
    """ Unpaired Space/Unpaired Time Peak Bias (%)"""
    return ((mod.max(axis = axis)-obs.max(axis = axis))/obs.max(axis = axis))*100.

def USUTPE(obs, mod, axis = None):
    """ Unpaired Space/Unpaired Time Peak Error (%)"""
    return (np.ma.abs(mod.max(axis = axis)-obs.max(axis = axis))/obs.max(axis = axis))*100.

def MNPB(obs, mod, paxis, axis = None):
    """ Mean Normalized Peak Bias (%)"""
    return ((mod.max(axis = paxis)-obs.max(axis = paxis))/obs.max(axis = paxis)).mean(axis = axis)*100.

def MdnNPB(obs, mod, paxis, axis = None):
    """ Median Normalized Peak Bias (%)"""
    return np.ma.median((mod.max(axis = paxis)-obs.max(axis = paxis))/obs.max(axis = paxis), axis = axis)*100.

def MNPE(obs, mod, paxis, axis = None):
    """ Mean Normalized Peak Error (%)"""
    return ((np.ma.abs(mod.max(axis = paxis)-obs.max(axis = paxis)))/obs.max(axis = paxis)).mean(axis = axis)*100.

def MdnNPE(obs, mod, paxis, axis = None):
    """ Median Normalized Peak Bias (%)"""
    return np.ma.median((np.ma.abs(mod.max(axis = paxis)-obs.max(axis = paxis)))/obs.max(axis = paxis), axis = axis)*100.

def NMPB(obs, mod, paxis, axis = None):
    """ Normalized Mean Peak Bias (%)"""
    return (mod.max(axis = paxis)-obs.max(axis = paxis)).mean(axis = axis)/obs.max(axis = paxis).mean(axis = axis)*100.

def NMdnPB(obs, mod, paxis, axis = None):
    """ Normalized Median Peak Bias (%)"""
    return np.ma.median((mod.max(axis = paxis)-obs.max(axis = paxis)), axis = axis)/np.ma.median(obs.max(axis = paxis), axis = axis)*100.

def NMPE(obs, mod, paxis, axis = None):
    """ Normalized Mean Peak Error (%)"""
    return (np.ma.abs(mod.max(axis = paxis)-obs.max(axis = paxis))).mean(axis = axis)/obs.max(axis = paxis).mean(axis = axis)*100.

def NMdnPE(obs, mod, paxis, axis = None):
    """ Normalized Median Peak Bias (%)"""
    return np.ma.median(np.ma.abs(mod.max(axis = paxis)-obs.max(axis = paxis)), axis = axis)/np.ma.median(obs.max(axis = paxis), axis = axis)*100.

def PSUTMNPB(obs, mod, axis = None):
    """ Paired Space/Unpaired Time Mean Normalized Peak Bias (%)"""
    return MNPB(obs, mod, paxis = 0, axis = None)

def PSUTMdnNPB(obs, mod, axis = None):
    """ Paired Space/Unpaired Time Median Normalized Peak Bias (%)"""
    return MdnNPB(obs, mod, paxis = 0, axis = None)

def PSUTMNPE(obs, mod, axis = None):
    """ Paired Space/Unpaired Time Mean Normalized Peak Error (%)"""
    return MNPE(obs, mod, paxis = 0, axis = None)

def PSUTMdnNPE(obs, mod, axis = None):
    """ Paired Space/Unpaired Time Median Normalized Peak Error (%)"""
    return MdnNPE(obs, mod, paxis = 0, axis = None)

def PSUTNMPB(obs, mod, axis = None):
    """ Paired Space/Unpaired Time Normalized Mean Peak Bias (%)"""
    return NMPB(obs, mod, paxis = 0, axis = None)

def PSUTNMPE(obs, mod, axis = None):
    """ Paired Space/Unpaired Time Normalized Mean Peak Error (%)"""
    return NMPE(obs, mod, paxis = 0, axis = None)

def PSUTNMdnPB(obs, mod, axis = None):
    """ Paired Space/Unpaired Time Normalized Median Peak Bias (%)"""
    return NMdnPB(obs, mod, paxis = 0, axis = None)

def PSUTNMdnPE(obs, mod, axis = None):
    """ Paired Space/Unpaired Time Normalized Median Peak Error (%)"""
    return NMdnPE(obs, mod, paxis = 0, axis = None)

def R2(obs, mod, axis = None):
    """ Coefficient of Determination (unit squared)"""
    from scipy.stats import pearsonr
    if axis is None:
        obsc, modc = matchedcompressed(obs, mod)
        return pearsonr(obsc, modc)[0]**2
    else:
        raise ValueError('Not ready yet')

def RMSE(obs, mod, axis = None):
    """ Root Mean Square Error (model unit)"""
    return np.ma.sqrt(((mod-obs)**2).mean(axis = axis))

def WDRMSE(obs, mod, axis = None):
    """ Wind Direction Root Mean Square Error (model unit)"""
    return np.ma.sqrt(((circlebias(mod-obs))**2).mean(axis = axis))

def RMSEs(obs, mod, axis = None):
    """Root Mean Squared Error (obs, mod_hat)"""
    if axis is None:
        try:
            from scipy.stats import linregress
            obsc, modc = matchedcompressed(obs, mod)
            m, b, rval, pval, stderr = linregress(obsc, modc)
            mod_hat = b + m * obs
            return RMSE(obs, mod_hat)
        except ValueError:
            return None
    else:
        raise ValueError('Not ready yet')

def matchmasks(a1, a2):
    mask = np.ma.getmaskarray(a1) | np.ma.getmaskarray(a2)
    return np.ma.masked_where(mask, a1), np.ma.masked_where(mask, a2)

def matchedcompressed(a1, a2):
    a1, a2 = matchmasks(a1, a2)
    return a1.compressed(), a2.compressed()

def RMSEu(obs, mod, axis = None):
    """Root Mean Squared Error (mod_hat, mod)"""
    if axis is None:
        try:
            from scipy.stats import linregress
            obsc, modc = matchedcompressed(obs, mod)
            m, b, rval, pval, stderr = linregress(obsc, modc)
            mod_hat = b + m * obs
            return RMSE(mod_hat, mod)
        except ValueError:
            return None
    else:
        raise ValueError('Not ready yet')

def d1(obs, mod, axis = None):
    """ Modified Index of Agreement, d1"""
    return 1.0 - (np.ma.abs(obs-mod)).sum(axis = axis)/(np.ma.abs(mod-obs.mean(axis = axis))+np.ma.abs(obs-obs.mean(axis = axis))).sum(axis = axis)

def E1(obs, mod, axis = None):
    """ Modified Coefficient of Efficiency, E1"""
    return 1.0 - (np.ma.abs(obs-mod)).sum(axis = axis)/(np.ma.abs(obs-obs.mean(axis = axis))).sum(axis = axis)

def IOA(obs, mod, axis = None):
    """ Index of Agreement, IOA"""
    obsmean = obs.mean(axis = axis)
    if not axis is None:
        obsmean = np.expand_dims(obsmean, axis = axis)
    return 1.0 - (np.ma.abs(obs-mod)**2).sum(axis = axis)/((np.ma.abs(mod-obsmean)+np.ma.abs(obs-obsmean))**2).sum(axis = axis)

def circlebias(b):
    b = np.ma.where(b > 180, b - 360, b)
    b = np.ma.where(b < -180, b + 360, b)
    return b

def WDIOA(obs, mod, axis = None):
    """ Wind Direction Index of Agreement, IOA"""
    obsmean = obs.mean(axis = axis)
    if not axis is None:
        obsmean = np.expand_dims(obsmean, axis = axis)
    b = circlebias(mod - obs)

    bhat = circlebias(mod - obsmean)

    ohat = circlebias(obs - obsmean)

    return 1.0 - (np.ma.abs(b)**2).sum(axis = axis)/((np.ma.abs(bhat)+np.ma.abs(ohat))**2).sum(axis = axis)

def AC(obs, mod, axis = None):
    """ Anomaly Correlation """
    obs_bar = obs.mean(axis = axis)
    if not axis is None:
        obs_bar = np.expand_dims(obs_bar, axis = axis)
    p1 = ((mod - obs_bar) * (obs - obs_bar)).sum(axis = axis)
    p2 = (((mod - obs_bar)**2).sum(axis = axis) * ((obs - obs_bar)**2).sum(axis = axis))**0.5
    return p1 / p2

def WDAC(obs, mod, axis = None):
    """ Wind Direction Anomaly Correlation """
    obs_bar = obs.mean(axis = axis)
    if not axis is None:
        obs_bar = np.expand_dims(obs_bar, axis = axis)
    p1 = (circlebias(mod - obs_bar) * circlebias(obs - obs_bar)).sum(axis = axis)
    p2 = ((circlebias(mod - obs_bar)**2).sum(axis = axis) * (circlebias(obs - obs_bar)**2).sum(axis = axis))**0.5
    return p1 / p2

def stat_spatial(ifile0, ifile1, funcs = __all__, variables = ['O3'], counties = False):
    """
    ifile0 - left hand side of equation
    ifile1 - right hand side of equation
    variables - list of variables to plot
    """
    from mpl_toolkits.basemap import Basemap
    from pylab import figure, show, colorbar
    lonlatcoords = getattr(ifile0, 'lonlatcoords', getattr(ifile1, 'lonlatcoords', ''))
    lon, lat = np.array(map(lambda x: map(float, x.split(',')), lonlatcoords.split('/'))).T
    latmin, latmax = lat.min(), lat.max()
    lonmin, lonmax = lon.min(), lon.max()
    bmap = Basemap(llcrnrlat = latmin, llcrnrlon = lonmin, urcrnrlat = latmax, urcrnrlon = lonmax, resolution = 'i')
    for vark in variables:
        for statname in funcs:
            statfunc = eval(statname)
            fig = figure()
            ax = fig.add_subplot(111)
            ax.set_title(vark + ' ' + statname)
            var_0 = ifile0.variables[vark]
            var_1 = ifile1.variables[vark]
            try:
                pidx = list(var_0.dimensions).index('points')
            except:
                pidx = list(var_1.dimensions).index('points')
            statvs = []
            for sitei in range(var_0.shape[pidx]):
                val_0 = var_0[:].take([sitei], axis = pidx).ravel()
                val_1 = var_1[:].take([sitei], axis = pidx).ravel()
                statvs.append(statfunc(val_0, val_1))
            dots = bmap.scatter(x = lon, y = lat, c = statvs, ax = ax, cmap = 'jet')
            
            bmap.drawcoastlines(ax = ax)
            bmap.drawcountries(ax = ax)
            bmap.drawstates(ax = ax)
            fig.colorbar(dots)
            if counties: bmap.counties(ax = ax)
            show()

def stat_timeseries(ifile0, ifile1, variables = ['O3'], counties = False):
    from pylab import figure, show, colorbar
    for vark in variables:
        for statname in __all__:
            statfunc = eval(statname)
            fig = figure()
            ax = fig.add_subplot(111)
            ax.set_title(vark + ' ' + statname)
            var_0 = ifile0.variables[vark]
            var_1 = ifile1.variables[vark]
            tidx = list(var_0.dimensions).index('time')
            statvs = []
            for timei in range(var_0.shape[tidx]):
                val_0 = var_0[:].take([timei], axis = tidx).ravel()
                val_1 = var_1[:].take([timei], axis = tidx).ravel()
                statvs.append(statfunc(val_0, val_1))
            dots = ax.plot(statvs)
            show() 

def main():
    from warnings import warn
    from pncparse import pncparse
    from PseudoNetCDF.core._functions import pncfunc, pncbfunc
    from PseudoNetCDF.pncparse import getparser, pncparse

    parser = getparser(has_ofile = False, plot_options = False, interactive = True)
    parser.add_argument('--funcs', default = __all__, type = lambda x: x.split(','), help='Functions to evaluate split by , (default: %s)' % ','.join(__all__))

    import numpy as np
    np.seterr(divide = 'ignore', invalid = 'ignore')
    ifiles, options = pncparse(has_ofile = False, interactive = True, parser = parser)
    if options.variables is None:
        options.variables = set(ifiles[0].variables.keys()).difference(options.coordkeys)
    console = createconsole(ifiles, options)
    warn("Assumes input order is obs model")
    ifile1, ifile2 = ifiles
    for k in options.funcs:
        console.locals[k] = func = eval(k)
        console.locals[k+'_f'] = ofile = pncbfunc(func, ifile1, ifile2)
        times = ifile1.variables['time']
        tstart = times[:].min()
        tstop = times[:].max()
        dt = tstop - tstart
        for vk in options.variables:
            if vk in ('time', 'TFLAG'): continue
            print '%.2f,%.2f,%.2f,%s,%s,%s,%f' % (tstart, tstop, dt, vk, func.__doc__.strip(), k, ofile.variables[vk].ravel()[0])
    np.seterr(divide = 'warn', invalid = 'warn')
    if options.interactive:
        console.interact()
if __name__ == '__main__':
    main()