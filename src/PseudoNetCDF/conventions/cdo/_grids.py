import textwrap
def griddes(infile):
    gridparams['GTYPE'] = 'unstructured'
    gridparams['NX'] = len(infile.dimensions['longitude'])
    gridparams['NY'] = len(infile.dimensions['latitude'])
    gridparams['XUNIT'] = 'degrees_east'
    gridparams['YUNIT'] = 'degrees_north'
    exec('NCELLS = NX * NY', None, gridparams)
    LONSTR = textwrap.wrap(' '.join(['%f' % lon for lon in infile.variables['longitude'][:]]), 80)
    LATSTR = textwrap.wrap(' '.join(['%f' % lat for lat in infile.variables['latitude'][:]]), 80)
    out = """
    gridtype  = %(GTYPE)s
    gridsize  = %(NCELLS)d
    xsize     = %(NX)d
    ysize     = %(NY)d
    xunits    = %(XUNIT)s
    yunits    = %(YUNIT)s
    xvals     = %(LONSTR)s
    yvals     = %(LATSTR)s
    """ % gridparams
    return out