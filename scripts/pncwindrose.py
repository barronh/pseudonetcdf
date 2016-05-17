#!/usr/bin/env python
# Run this script with pnc options
import numpy as np
def pncwindrose(ifiles, args):
    if len(ifiles) > 2:
        raise ValueError('Requires one or two input files (received %d). To combine, use --stack' % len(ifiles))
        
    
    ws = np.ma.masked_invalid(eval(args.windspeed, None, ifiles[0].variables)[:])
    wd = np.ma.masked_invalid(eval(args.winddir, None, ifiles[0].variables)[:])
    if len(ifiles) == 2:
        maskvar = eval(args.maskvar, None, ifiles[1].variables)
        ws = np.ma.masked_where(maskvar, ws[:])
        wd = np.ma.masked_where(maskvar, wd[:])
        
    windrose(wd[:], ws[:], args)

def windrose(wd, ws, args):
    from matplotlib import pyplot as plt
    plt.rcParams['text.usetex'] = False
    width = args.binwidth
    halfwidth = width / 2.
    fig = plt.figure(figsize = (6, 7))
    ax = fig.add_axes([0.1, 0.2, 0.8, 0.7], polar = True)
    ax.set_theta_offset(np.radians(90))
    ax.set_theta_direction(-1)
    ubs = args.bounds[1:]
    lbs = args.bounds[:-1]
    ubcs = plt.cm.jet(np.arange(len(ubs), dtype = 'f') / ubs.size)
    color_dict = dict(zip(ubs, ubcs))
    mask = np.ma.getmaskarray(ws) | np.ma.getmaskarray(wd) | (ws < args.bounds[0]) | (ws > args.bounds[-1])
    ws = np.ma.masked_where(mask, ws).compressed()
    wd = np.ma.masked_where(mask, wd).compressed()
    fig.suptitle('%s Wind Rose (%d)' % (args.title,ws.size), size = 20)
    if args.fromnorth:
        # NNE, WNE, WSE, SSE, SSW, WSW, WNW, NNW
        thetas = ((np.radians((wd) // width * width)).astype('f') % (2 * np.pi)).astype('f')
    else:
        # N, NE, E, SE, S, SW, W, NW
        thetas = (np.radians((wd.astype('d') + halfwidth) // width * width - halfwidth) % (2 * np.pi)).astype('f')
    bounds = np.degrees(thetas)[:, None] + np.array([0, width])[None, :] % 360
    uthetas = np.unique(thetas)
    nwinds = float(ws.size)
    for theta in uthetas:
        tws = ws[thetas == theta]
        for ubi, ub in enumerate(ubs[::-1]):
            ubc = color_dict[ub]
            pct = (tws < ub).sum() / nwinds * 100.
            ax.bar(theta, pct, width = np.radians(width), bottom = 0.0, color = ubc)
    labels = [(('<%s m/s' % ub) if ub != np.inf else ('>%s m/s' % lb)) for lb, ub in zip(lbs, ubs)]
    rectangles = [plt.Rectangle((0, 0), 1, 1, color = ubc) for ubc in ubcs]
    fig.legend(rectangles, labels, ncol = 3, bbox_to_anchor = (.5, 0.025), loc = 'lower center')
    ticks = np.linspace(0, args.maxpct, 6)[1:]
    labels = ['%s%%' % i for i in ticks]
    plt.yticks(ticks, labels)
    ax.set_rmax(args.maxpct)
    ax.set_clip_on(False)
    figpath = args.outpath
    fig.savefig(figpath)
    if args.verbose > 0: print('Saved fig', figpath)
    plt.close(fig)
    
if __name__ == '__main__':
    from PseudoNetCDF.pncparse import getparser, pncparse
    parser = getparser(has_ofile = True, plot_options = False, interactive = False)
    parser.add_argument('--wind-speed', dest = 'windspeed', default = 'WS', help = 'Wind speed variable')
    parser.add_argument('--mask-var', dest = 'maskvar', default = 'WS', help = 'Variable from file two for masking')
    parser.add_argument('--wind-dir', dest = 'winddir', default = 'WD', help = 'Wind direction variable')
    parser.add_argument('--title', default = None, help = 'Title for wind rose')
    parser.add_argument('--binwidth', dest = 'binwidth', type = float, default = 45., help = 'Direction bin widths in degrees.')
    parser.add_argument('--max-pct', dest = 'maxpct', type = float, default = 50., help = 'Maximum percent on plot')
    parser.add_argument('--fromnorth', dest = 'fromnorth', action = 'store_true', default = False, help = 'Bins should start from north (i.e., not bound north)')
    parser.add_argument('--bounds', dest = 'bounds', type = lambda x: np.array(eval(x)), default = np.array([1, 2, 3, 5, 10, np.inf]), help = 'Boundaries for wind speed')
    parser.epilog = """
Example:
    $ pncwindrose.py -s time,600,None --max-pct=35 --wind-dir wd --wind-speed ws  --mask-var="ws[:].mask.all(0)[None,:].repeat(ws.shape[0], 0)" --bounds="[0.5,1,2,3.5,4.5,6.5,np.inf]"  --title obs  Bogota_2012_s12345.nc df.nc oct_obs

Description:
    Create windrose for timestep 600 and on where wind-direction is name wd, and windspeed is named ws. Masking is applied from df.nc file where all times for ws are masked (useful for ignoring monitors). Values less than 0.5 will be masked.
"""
    ifiles, args = pncparse(has_ofile = True, plot_options = False, interactive = False, args = None, parser = parser)
    out = pncwindrose(ifiles, args)
    