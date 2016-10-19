#!/usr/bin/env python
def main():
    from PseudoNetCDF.pncparse import getparser, pncparse
    from PseudoNetCDF.plotutil.pncmap import makemaps
    parser = getparser(has_ofile = True, map_options = True, plot_options = True, interactive = True)
    ifiles, args = pncparse(has_ofile = True, plot_options = True, interactive = True, parser = parser)
    args.ifiles = ifiles
    makemaps(args)

if __name__ == '__main__':
    main()
