#!/usr/bin/env python
from PseudoNetCDF.plotutil.vertprofile import vertprofileplot, add_vertprofile_options
from PseudoNetCDF.pncparse import pncparse, getparser

if __name__ == '__main__':
    vertparser = getparser(has_ofile = True, plot_options = True, interactive = False)
    add_vertprofile_options(vertparser)
    ifiles, args = pncparse(has_ofile = True, parser = vertparser)
    vertprofileplot(ifiles, args)
