__all__ = ['register', 'get_cmap', 'register_norm', 'get_norm']

from matplotlib.colors import from_levels_and_colors, Normalize
from matplotlib.cm import get_cmap
import matplotlib as mpl

try:
    register = mpl.colormaps.register
except AttributeError:
    register = mpl.cm.register_cmap


_registered_norms = {}


def get_norm(name):
    global _registered_norms
    return _registered_norms[name]


def register_norm(name, norm):
    global _registered_norms
    if isinstance(norm, Normalize):
        _registered_norms[name] = norm
    else:
        raise TypeError(
            'cmap should be an instance of Normalize; got ' + str(type(norm)))


def _addusepaaqi():
    aqicolors = ['green', 'yellow', 'orange', 'red', 'purple', 'maroon']
    aqic = aqicolors
    bounds_colors = {}
    bounds_colors[('O3', 'mda8', 'ppbv')] = ([0, 55, 71, 86, 106, 200],
                                             aqic[:5] + ['black'], 'max')
    bounds_colors[('O3', 'mda8', 'ppmv')] = ([0, 0.055, 0.071, 0.086,
                                              0.106, 0.0200],
                                             aqic[:5] + ['black'], 'max')
    bounds_colors[('O3', 'mda1', 'ppbv')] = ([125, 165, 205, 405, 604],
                                             ['white'] + aqic[2:] + ['black'],
                                             'both')
    pmb = [0, 12.1, 35.5, 55.5, 150.5, 250.5, 500.4]
    pmc = aqic + ['black']
    bounds_colors[('PM25', 'a24', 'micrograms/m**3')] = (pmb, pmc, 'max')
    for key, (bnds, colors, extend) in bounds_colors.items():
        print(key)
        name = '_'.join(('usepa', 'aqi',) + key)
        cmap, norm = from_levels_and_colors(bnds, colors, extend=extend)
        register(name=name, cmap=cmap)
        register_norm(name=name, norm=norm)


_addusepaaqi()


if __name__ == '__main__':
    from matplotlib import pyplot as plt
    import numpy as np
    varkey = 'O3'
    intv = 'mda8'
    unit = 'ppbv'

    colorscale = 'usepa_aqi_' + '_'.join([varkey, intv, unit])
    norm = get_norm(colorscale)
    cmap = get_cmap(colorscale)
    data = np.linspace(25, 200, 16).reshape(4, 4)
    print(np.flipud(data))
    plt.pcolor(data, norm=norm, cmap=cmap)
    plt.savefig('test.png')
