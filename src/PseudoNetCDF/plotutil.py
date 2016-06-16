import pylab
pl = pylab
import matplotlib.pyplot as pyplot
plt = pyplot

Normalize = pl.matplotlib.colors.Normalize
LogNorm = pl.matplotlib.colors.LogNorm
BoundaryNorm = pl.matplotlib.colors.BoundaryNorm

def SegmentedNorm(vmin, vmax, bins = 10, ncolors = 256):
    boundaries = np.linspace(vmin, vmax, bins + 1)
    return BoundaryNorm(boundaries, ncolors)

def SegmentedLogNorm(vmin, vmax, bins = 10, ncolors = 256, full_levels = False):
    """
    Create BoundaryNorm with (N=bins) log spacing.  If full_levels,
    then colorbar will start on a full log level and end on a full
    log level.
    
    
    
    boundaries = np.logspace(np.log10(vmin), np.log10(vmax), bins + 1)
    return BoundaryNorm(boundaries, ncolors)
    """
    if full_levels:
        boundaries = np.logspace(np.floor(np.log10(vmin)), np.ceil(np.log10(vmax)), bins + 1)
    else:
        boundaries = np.logspace(np.log10(vmin), np.log10(vmax), bins + 1)
    return BoundaryNorm(boundaries, ncolors)

LogFormatter = pl.matplotlib.ticker.LogFormatter
ScalarFormatter = pl.matplotlib.ticker.ScalarFormatter
