__all__ = ['NO', 'NP', 'NOP', 'MO', 'MdnO', 'MP', 'MdnP', 'RM', 'RMdn', 'MB', 'MdnB', 'FB', 'MNB', 'MdnNB', 'NMB', 'NMdnB', 'USUTPB', 'PSUTMNPB', 'PSUTMdnNPB', 'PSUTNMPB', 'PSUTNMdnPB', 'ME', 'MdnE', 'FE', 'MNE', 'MdnNE', 'NME', 'NMdnE', 'USUTPE',  'PSUTMNPE', 'PSUTMdnNPE', 'PSUTNMPE', 'PSUTNMdnPE', 'R2', 'RMSE', 'RMSEs', 'RMSEu', 'E1', 'IOA', 'd1']
from pncload import createconsole
import numpy as np
def MNB(obs, mod):
    """ Mean Normalized Bias (%)"""
    return  np.ma.masked_invalid((mod-obs)/obs).mean()*100.

def MNE(obs, mod):
    """ Mean Normalized Gross Error (%)"""
    return np.ma.masked_invalid(np.ma.abs(mod-obs)/obs).mean()*100.

def MdnNB(obs, mod):
    """ Median Normalized Bias (%)"""
    return np.ma.median(np.ma.masked_invalid((mod-obs)/obs))*100.

def MdnNE(obs, mod):
    """ Median Normalized Gross Error (%)"""
    return np.ma.median(np.ma.masked_invalid(np.ma.abs(mod-obs)/obs))*100.

def NMdnE(obs, mod):
    """ Normalized Median Gross Error (%)"""
    return np.ma.masked_invalid(np.ma.abs(mod-obs).mean()/obs.mean())*100.

def NO(obs, mod):
    """ N Observations (#)"""
    return (np.ma.getmaskarray(obs) == False).sum()

def NOP(obs, mod):
    """ N Observations/Prediction Pairs (#)"""
    obsc, modc = matchedcompressed(obs, mod)
    return obsc.size

def NP(obs, mod):
    """ N Predictions (#)"""
    return (np.ma.getmaskarray(mod) == False).sum()


def MO(obs, mod):
    """ Mean Observations (obs unit)"""
    return obs.mean()

def MP(obs, mod):
    """ Mean Predictions (model unit)"""
    return mod.mean()

def MdnO(obs, mod):
    """ Median Observations (obs unit)"""
    return np.ma.median(obs)

def MdnP(obs, mod):
    """ Median Predictions (model unit)"""
    return np.ma.median(mod)

def RM(obs, mod):
    """ Mean Ratio Observations/Predictions (none)"""
    return np.ma.masked_invalid(obs/mod).mean()

def RMdn(obs, mod):
    """ Median Ratio Observations/Predictions (none)"""
    return np.ma.median(np.ma.masked_invalid(obs/mod))

def MB(obs, mod):
    """ Mean Bias"""
    return (mod-obs).mean()

def MdnB(obs, mod):
    """ Median Bias"""
    return np.ma.median(mod-obs)

def NMB(obs, mod):
    """ Normalized Mean Bias (%)"""
    return (mod-obs).sum()/obs.sum()*100.

def NMdnB(obs, mod):
    """ Normalized Median Bias (%)"""
    return np.ma.median(mod-obs)/np.ma.median(obs)*100.

def FB(obs, mod):
    """ Fractional Bias (%)"""
    obsc, modc = matchedcompressed(obs, mod)
    return ((np.ma.masked_invalid((mod-obs)/(mod+obs))).mean()*2.)*100.

def ME(obs, mod):
    """ Mean Gross Error (model and obs unit)"""
    return np.ma.abs(mod-obs).mean()

def MdnE(obs, mod):
    """ Median Gross Error (model and obs unit)"""
    return np.ma.median(np.ma.abs(mod-obs))

def NME(obs, mod):
    """ Normalized Mean Error (%)"""
    out = (np.ma.abs(mod-obs).sum()/obs.sum())*100
    return out

def NMdnE(obs, mod):
    """ Normalized Median Error (%)"""
    out = np.ma.median(np.ma.abs(mod-obs))/np.ma.median(obs)*100
    return out

def FE(obs, mod):
    """ Fractional Error (%)"""
    obsc, modc = matchedcompressed(obs, mod)
    return (np.ma.abs(mod-obs)/(mod+obs)).mean()*2.*100.
    #Part 1b. For daily peak (or some period)

def USUTPB(obs, mod):
    """ Unpaired Space/Unpaired Time Peak Bias (%)"""
    return ((mod.max()-obs.max())/obs.max())*100.

def USUTPE(obs, mod):
    """ Unpaired Space/Unpaired Time Peak Error (%)"""
    return (np.ma.abs(mod.max()-obs.max())/obs.max())*100.

def PSUTMNPB(obs, mod):
    """ Paired Space/Unpaired Time Mean Normalized Peak Bias (%)"""
    return ((mod.max(0)-obs.max(0))/obs.max(0)).mean()*100.

def PSUTMdnNPB(obs, mod):
    """ Paired Space/Unpaired Time Median Normalized Peak Bias (%)"""
    return np.ma.median((mod.max(0)-obs.max(0))/obs.max(0))*100.

def PSUTMNPE(obs, mod):
    """ Paired Space/Unpaired Time Mean Normalized Peak Error (%)"""
    return (np.ma.abs(mod.max(0)-obs.max(0))/obs.max(0)).mean()*100.

def PSUTMdnNPE(obs, mod):
    """ Paired Space/Unpaired Time Median Normalized Peak Error (%)"""
    return np.ma.median(np.ma.abs(mod.max(0)-obs.max(0))/obs.max(0))*100.

def PSUTNMPB(obs, mod):
    """ Paired Space/Unpaired Time Normalized Mean Peak Bias (%)"""
    return (mod.max(0)-obs.max(0)).mean()/obs.max(0).mean()*100.

def PSUTNMPE(obs, mod):
    """ Paired Space/Unpaired Time Normalized Mean Peak Error (%)"""
    return np.ma.abs(mod.max(0)-obs.max(0)).mean()/obs.max(0).mean()*100.

def PSUTNMdnPB(obs, mod):
    """ Paired Space/Unpaired Time Normalized Median Peak Bias (%)"""
    return np.ma.median(mod.max(0)-obs.max(0))/np.ma.median(obs.max(0))*100.

def PSUTNMdnPE(obs, mod):
    """ Paired Space/Unpaired Time Normalized Median Peak Error (%)"""
    return np.ma.median(np.ma.abs(mod.max(0)-obs.max(0)))/np.ma.median(obs.max(0))*100.

def R2(obs, mod):
    """ Correlation Coefficient (unit squared)"""
    from scipy.stats import pearsonr
    obsc, modc = matchedcompressed(obs, mod)
    return pearsonr(obsc, modc)[0]**2

def RMSE(obs, mod):
    """ Root Mean Square Error (model unit)"""
    return np.ma.sqrt(((mod-obs)**2).mean())

def RMSEs(obs, mod):
    """Root Mean Squared Error (obs, mod_hat)"""
    try:
        from scipy.stats import linregress
        obsc, modc = matchedcompressed(obs, mod)
        m, b, rval, pval, stderr = linregress(obsc, modc)
        mod_hat = b + m * obs
        return RMSE(obs, mod_hat)
    except ValueError:
        return None

def matchmasks(a1, a2):
    mask = np.ma.getmaskarray(a1) | np.ma.getmaskarray(a2)
    return np.ma.masked_where(mask, a1), np.ma.masked_where(mask, a2)

def matchedcompressed(a1, a2):
    a1, a2 = matchmasks(a1, a2)
    return a1.compressed(), a2.compressed()

def RMSEu(obs, mod):
    """Root Mean Squared Error (mod_hat, mod)"""
    try:
        from scipy.stats import linregress
        obsc, modc = matchedcompressed(obs, mod)
        m, b, rval, pval, stderr = linregress(obsc, modc)
        mod_hat = b + m * obs
        return RMSE(mod_hat, mod)
    except ValueError:
        return None

def d1(obs, mod):
    """ Modified Index of Agreement, d1"""
    return 1.0 - (np.ma.abs(obs-mod)).sum()/(np.ma.abs(mod-obs.mean())+np.ma.abs(obs-obs.mean())).sum()

def E1(obs, mod):
    """ Modified Coefficient of Efficiency, E1"""
    return 1.0 - (np.ma.abs(obs-mod)).sum()/(np.ma.abs(obs-obs.mean())).sum()

def IOA(obs, mod):
    """ Index of Agreement, IOA"""
    return 1.0 - (np.ma.abs(obs-mod)**2).sum()/((np.ma.abs(mod-obs.mean())+np.ma.abs(obs-obs.mean()))**2).sum()

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
    from pncparse import pncparser
    from PseudoNetCDF.core._functions import pncfunc, pncbfunc
    import numpy as np
    np.seterr(divide = 'ignore', invalid = 'ignore')
    ifiles, options = pncparser(has_ofile = False)
    console = createconsole(ifiles, options)
    warn("Assumes input order is obs model")
    ifile1, ifile2 = ifiles
    for k in __all__:
        console.locals[k] = func = eval(k)
        console.locals[k+'_f'] = ofile = pncbfunc(func, ifile1, ifile2)
        for vk in ofile.variables.keys():
            print vk, func.__doc__.strip(), '(' + k + ')', ofile.variables[vk].ravel()[0]
    np.seterr(divide = 'warn', invalid = 'warn')
    if options.interactive:
        console.interact()
if __name__ == '__main__':
    main()