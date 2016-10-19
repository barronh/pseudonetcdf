#!/usr/bin/env python

if __name__ == '__main__':
    from PseudoNetCDF.plotutil.pncscatter import pncscatter
    from PseudoNetCDF.pncparse import getparser, pncparse
    parser = getparser(plot_options = True, has_ofile = True)
    parser.epilog += """
    -----
pncscatter inobs inmod target [target ...]
inobs - path to obs (or other model)
inmod - path to mod (or other obs)
target - variable name
"""
    ifiles, args = pncparse(plot_options = True, has_ofile = True, parser = parser)
    pncscatter(args)
