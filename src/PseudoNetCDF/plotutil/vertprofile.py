#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from warnings import warn
from netCDF4 import Dataset
from collections import defaultdict
from datetime import datetime, timedelta
from glob import glob

from PseudoNetCDF.pncparse import AggCommaString

unitconvert = {('ppmV', 'ppb'): lambda x: x * 1000.}


def add_vertprofile_options(vertparser):
    vpadd = vertparser.add_argument
    vpadd("--sigma", dest="sigma", action="store_true", default=False,
          help="Plot on sigma coordinate instead of pressure")

    vpadd("--scale", dest="scale", type=str, default='log',
          help="Defaults to log, but linear and semilog are also options.")

    vpadd("--minmax", dest="minmax", type=str, default="None,None",
          help="Use values to set range (xmin, xmax); defaults None,None.")

    vpadd("--mask-zeros", dest="maskzeros", action="store_true",
          default=False, help="Defaults False.")

    vpadd("--minmaxq", dest="minmaxq", type=str, default='0,100',
          help="Use quartiles to set range (xmin, xmax); defaults 0,100.")

    vpadd("--out-unit", dest="outunit", type=str, default='ppb',
          help="Defaults ppb.")

    vpadd("--tes-paths", dest="tespaths", type=str,
          action=AggCommaString, default='',
          help="Plot tes on top of boundary from paths; defaults to []")

    vpadd("--omi-paths", dest="omipaths", type=str,
          default='', action=AggCommaString,
          help="Plot omi on top of boundary from paths; defaults to []")

    vpadd("--itertime", dest="itertime", default=False, action='store_true',
          help="Iterate over times and plot each one.")

    vpadd("--edges", dest="edges", default=False, action="store_true",
          help="Plot S,E,N,W edges instead of a single plot.")

# vertparser = getparser(has_ofile = True, plot_options = True,
#                        interactive = False)
# add_vertprofile_options(vertparser)
# vertparser.epilog = """
# Example:
# $ pncvertprofile.py outputs/ts20120301.bpch.BCON.nc test_profile -v O3 \
#    --edges #--minmaxq .5,99.5 --tes-paths=~/Data/test/*
# """


def maskfilled(v):
    return np.ma.masked_values(v[:], v.attrs['_FillValue'])


def plot_omi(ax, lon_bnds, lat_bnds, omipaths, key='O3Profile',
             airden=None, airdenvert=None):
    import h5py
    # from scipy.constants import Avogadro
    allx = []
    ally = []
    lon_bnds = lon_bnds[:]
    lat_bnds = lat_bnds[:]
    if len(omipaths) == 0:
        return None, None
    outomipaths = []
    for omipath in omipaths:
        outomipaths.extend(glob(omipath))
    omipaths = outomipaths
    omipaths.sort()
    for path in omipaths:
        print(path)
        f = h5py.File(path, mode='r')
        swaths = f['HDFEOS']['SWATHS'][key]
        latsv = swaths['Geolocation Fields']['Latitude']
        lonsv = swaths['Geolocation Fields']['Longitude']
        lats = maskfilled(latsv).reshape(-1, 1)
        lons = maskfilled(lonsv).reshape(-1, 1)
        inboth = matchspace(lons, lats, lon_bnds, lat_bnds)
        if inboth.filled(False).sum(0) > 0:
            print('******** FOUND ******', path)
            pressurev = swaths['Geolocation Fields']['Pressure']
            altitudev = swaths['Geolocation Fields']['Altitude']
            nk = pressurev.shape[-1]
            speciesv = swaths['Data Fields']['O3']
            species = maskfilled(speciesv).reshape(-1, nk - 1)

            # Pressure is from top to bottom
            pressure = maskfilled(pressurev).reshape(-1, nk)
            altitude = maskfilled(altitudev).reshape(-1, nk)
            dz = -(altitude[:, 1:] - altitude[:, :-1]) * 1000. * 100.
            # Lacis 1990 1 DU = 2.69e16 molecules cm-2
            species = species * .001 / dz  # cm /cm # vrm
            x = pressure[inboth]
            y = species[inboth]
            allx.append(x)
            ally.append(y)
        else:
            wtmp = 'No data found for %s: lons (%s, %s) and lats (%s, %s); '
            wtmp += 'lonbs (%s, %s) and latbs (%s, %s);'
            wvals = (path, lons.min(), lons.max(), lats.min(), lats.max())
            wvals = wvals + (lon_bnds.min(), lon_bnds.max(), lat_bnds.min())
            wvals = wvals + (lat_bnds.max(),)
            warn(wtmp % wvals)

    if len(allx) == 0:
        print('*' * 80 + '\n\nNo OMI DATA FOUND AT ALL\n\n' + '*' * 80)
        return None, None

    var = np.ma.masked_values(np.ma.concatenate(ally, axis=0), -999.) * 1e9
    var = var.reshape(-1, var.shape[-1])
    vertcrd = np.ma.masked_values(np.ma.concatenate(
        allx, axis=0), -999.).mean(0)  # .reshape(-1, 2).mean(1)
    omil, omir = minmaxmean(ax, var.T.repeat(2, 0), vertcrd.repeat(2, 0)[1:-1],
                            ls='-', lw=2, color='b', facecolor='b',
                            edgecolor='b', alpha=.2, zorder=5, label='OMI')
    ax.text(.05, .7, 'OMI = %d' % var.shape[0], transform=ax.transAxes)
    return omil, omir


def matchspace(lons, lats, lon_bnds, lat_bnds):
    slicer = [slice(None, None)] + [None] * (lon_bnds.ndim - 1)
    inlon = np.logical_and(lons[slicer] >= lon_bnds[None, ..., 0],
                           lons[slicer] <= lon_bnds[None, ..., 1])
    inlat = np.logical_and(lats[slicer] >= lat_bnds[None, ..., 0],
                           lats[slicer] <= lat_bnds[None, ..., 1])
    inboth = np.logical_and(inlon, inlat).reshape(inlon.shape[0], -1).any(1)
    return inboth


def plot_tes(ax, lon_bnds, lat_bnds, tespaths):
    allx = []
    ally = []
    lon_bnds = lon_bnds[:]
    lat_bnds = lat_bnds[:]
    if len(tespaths) == 0:
        return None, None
    outtespaths = []
    for tespath in tespaths:
        outtespaths.extend(glob(tespath))
    tespaths = outtespaths
    tespaths.sort()
    for path in tespaths:
        f = Dataset(path)
        lats = f.variables['latitude'][:][:, None]
        lons = f.variables['longitude'][:][:, None]
        pressure = f.variables['pressure']
        species = f.variables['species']
        inboth = matchspace(lons, lats, lon_bnds, lat_bnds)
        if inboth.sum(0) > 0:
            print('******** FOUND ******', path)
            x = pressure[inboth]
            y = species[inboth]
            allx.append(x)
            ally.append(y)
        else:
            warn('No data found for %s' % path)

    if len(allx) == 0:
        return None, None
    var = np.ma.masked_values(np.ma.concatenate(ally, axis=0), -999.) * 1e9
    var = var.reshape(-1, var.shape[-1])
    vertcrd = np.ma.masked_values(
        np.ma.concatenate(allx, axis=0), -999.).mean(0)
    tesl, tesr = minmaxmean(ax, var.T, vertcrd, ls='-', lw=2, color='r',
                            facecolor='r', edgecolor='r', alpha=.2, zorder=2,
                            label='TES (%d)' % var.shape[0])
    return tesl, tesr


def minmaxmean(ax, vals, vertcrd, **kwds):
    minval = vals.min(1)
    meanval = vals.mean(1)
    maxval = vals.max(1)
    linekwds = kwds.copy()
    linekwds['color'] = linekwds.pop('facecolor')
    linekwds.pop('edgecolor')
    linekwds.pop('alpha')
    fillkwds = kwds.copy()
    fillkwds['ls'] = 'solid'

    line, = ax.plot(meanval, vertcrd, **linekwds)

    x = np.ma.concatenate([minval[:vertcrd.size], maxval[:vertcrd.size][::-1]])
    y = np.ma.concatenate([vertcrd[:], vertcrd[::-1]])
    mask = x.mask | y.mask
    x = np.ma.masked_where(mask, x).compressed()
    y = np.ma.masked_where(mask, y).compressed()
    range, = ax.fill(x, y, **fillkwds)
    return line, range


def plotprofile(args):
    vertprofileplot(args.ifiles, args)


def vertprofileplot(ifiles, args):
    if args.variables is None:
        raise ValueError('User must specify variable(s) to plot:\n%s' %
                         '\n\t'.join(ifiles[0].variables.keys()))
    from PseudoNetCDF.coordutil import getsigmamid, getpresmid, gettimes
    import matplotlib.pyplot as plt
    plt.rcParams['text.usetex'] = False
    # from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm
    # from matplotlib.colors import LogNorm
    scale = args.scale
    minmax = eval(args.minmax)
    minmaxq = eval(args.minmaxq)
    sigma = args.sigma
    maskzeros = args.maskzeros
    outunit = args.outunit
    tespaths = args.tespaths
    omipaths = args.omipaths
    edges = args.edges
    try:
        f, = ifiles
    except Exception:
        raise ValueError('curtain plot expects one file when done. ' +
                         'Try stack time --stack=time to concatenate')

    # Add CF conventions if necessary
    if 'latitude_bounds' not in f.variables.keys():
        try:
            from PseudoNetCDF import getvarpnc
            from PseudoNetCDF.conventions.ioapi import add_cf_from_ioapi
            f = getvarpnc(f, None)
            add_cf_from_ioapi(f)
        except Exception:
            pass
    if sigma:
        vertcrd = getsigmamid(f)
    else:
        vertcrd = getpresmid(f, pref=101325., ptop=getattr(f, 'VGTOP', 10000))
        if vertcrd.max() > 2000:
            vertcrd /= 100.

    try:
        lonb = f.variables['geos_longitude_bounds']
        latb = f.variables['geos_latitude_bounds']
    except Exception:
        lonb = f.variables['longitude_bounds']
        latb = f.variables['latitude_bounds']
    for var_name in args.variables:
        temp = defaultdict(lambda: 1)
        try:
            eval(var_name, None, temp)
            var = eval(var_name, None, f.variables)[:]
        except Exception:
            temp[var_name]
            var = f.variables[var_name][:]
        if maskzeros:
            var = np.ma.masked_values(var, 0)
        vkeys = [k for k in temp.keys()]
        unit = f.variables[vkeys[0]].units.strip()
        if unit in unitconvert:
            var = unitconvert.get((unit, outunit), lambda x: x)(var)
        else:
            outunit = unit
        vmin, vmax = np.percentile(
            np.ma.compressed(var).ravel(), list(minmaxq))
        if minmax[0] is not None:
            vmin = minmax[0]
        if minmax[1] is not None:
            vmax = minmax[1]
        if edges:
            fig = plt.figure(figsize=(16, 4))
            offset = 0.05
            ax = fig.add_axes([.1 - offset, .15, .22, .725])
            ax = fig.add_axes([.325 - offset, .15, .22, .725])
            ax = fig.add_axes([.55 - offset, .15, .22, .725])
            ax = fig.add_axes([.775 - offset, .15, .22, .725])
            ss = 0
            se = ss + f.NCOLS + 1
            es = se
            ee = se + f.NROWS + 1
            ns = ee
            ne = ee + f.NCOLS + 1
            ws = ne
            we = ws + f.NROWS + 1
            axs = fig.axes
            for ax in fig.axes[1:]:
                ax.yaxis.set_major_formatter(plt.NullFormatter())

            vars = [var[:, :, ss:se], var[:, :, es:ee], var[:, :, ns:ne]
                    [:, :, ::-1], var[:, :, ws:we][:, :, ::-1]]
            lonbss = [lonb[ss:se], lonb[es:ee],
                      lonb[ns:ne][::-1], lonb[ws:we][::-1]]
            latbss = [latb[ss:se], latb[es:ee],
                      latb[ns:ne][::-1], latb[ws:we][::-1]]

        else:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            axs = fig.axes
            vars = [var]
            ismultidimcoord = (lonb.dimensions == ('longitude', 'nv') and
                               latb.dimensions == ('latitude', 'nv'))
            if ismultidimcoord:
                lonbss = [lonb[:][None, :, :]]
                latbss = [latb[:][:, None, :]]
            else:
                lonbss = [lonb[:]]
                latbss = [latb[:]]
        for ax, var, lonbs, latbs in zip(axs, vars, lonbss, latbss):
            vals = var.swapaxes(0, 1).reshape(var.shape[1], -1)
            modl, modr = minmaxmean(ax, vals, vertcrd, facecolor='k',
                                    edgecolor='k', alpha=.2, zorder=4,
                                    label='mod (%d)' % vals.shape[1],
                                    ls='-', lw=2, color='k')
            llines = [(modl, modr)]
            ymin, ymax = vertcrd.min(), vertcrd.max()
            ax.set_ylim(ymax, ymin)
            ax.set_xscale(scale)
            ax.set_xlim(vmin, vmax)
            # if scale == 'log':
            #   ax.set_xticklabels(['%.1f' % (10**x) for x in ax.get_xticks()])

            if 'TFLAG' in f.variables.keys():
                SDATE = f.variables['TFLAG'][:][0, 0, 0]
                EDATE = f.variables['TFLAG'][:][-1, 0, 0]
                STIME = f.variables['TFLAG'][:][0, 0, 1]
                ETIME = f.variables['TFLAG'][:][-1, 0, 1]
                if SDATE == 0:
                    SDATE = 1900001
                    EDATE = 1900001
                sdate = datetime.strptime(
                    '%07d %06d' % (SDATE, STIME), '%Y%j %H%M%S')
                edate = datetime.strptime(
                    '%07d %06d' % (EDATE, ETIME), '%Y%j %H%M%S')
            elif 'tau0' in f.variables.keys():
                sdate = datetime(1985, 1, 1, 0) + \
                    timedelta(hours=f.variables['tau0'][0])
                edate = datetime(1985, 1, 1, 0) + \
                    timedelta(hours=f.variables['tau1'][-1])
            else:
                times = gettimes(f)
                sdate = times[0]
                edate = times[-1]

            if len(tespaths) > 0:
                tesl, tesr = plot_tes(ax, lonbs, latbs, tespaths)
                if tesl is not None:
                    llines.append((tesl, tesr))
            if len(omipaths) > 0:
                airdenvals = f.variables['AIRDEN'][:].mean(0).mean(1)
                omil, omir = plot_omi(ax, lonbs, latbs, omipaths,
                                      airden=airdenvals, airdenvert=vertcrd)
                if omil is not None:
                    llines.append((omil, omir))

        try:
            title = '%s to %s' % (sdate.strftime(
                '%Y-%m-%d'), edate.strftime('%Y-%m-%d'))
        except Exception:
            title = var_name
        if sigma:
            axs[0].set_ylabel('sigma')
        else:
            axs[0].set_ylabel('pressure')

        xmax = -np.inf
        xmin = np.inf
        for ax in fig.axes:
            tmp_xmin, tmp_xmax = ax.get_xlim()
            xmax = max(tmp_xmax, xmax)
            xmin = min(tmp_xmin, xmin)
        for ax in fig.axes:
            ax.set_xlim(xmin, xmax)

        if len(axs) == 1:
            axs[0].set_xlabel('%s %s' % (var_name, outunit))
        else:
            axs[0].set_xlabel('South')
            axs[1].set_xlabel('East')
            axs[2].set_xlabel('North')
            axs[3].set_xlabel('West')
            fig.text(.5, .90, '%s %s' % (var_name, outunit),
                     horizontalalignment='center', fontsize=16)
        nl = 0
        for ax in axs:
            if len(ax.get_lines()) > nl:
                nl = len(ax.get_lines())
                plt.sca(ax)

        llabels = [_l[0].get_label() for _l in llines]
        plt.legend(llines, llabels, bbox_to_anchor=(.1, 1),
                   loc='upper left', bbox_transform=fig.transFigure, ncol=6)
        if edges:
            fig.text(0.95, 0.975, title, horizontalalignment='right',
                     verticalalignment="top", fontsize=16)
        else:
            fig.text(0.95, 0.025, title, horizontalalignment='right',
                     verticalalignment="bottom", fontsize=16)
        figpath = args.outpath + var_name + '.' + args.figformat
        fig.savefig(figpath)
        if args.verbose > 0:
            print('Saved fig', figpath)
        # plt.close(fig)
    return fig
