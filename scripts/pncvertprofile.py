#!/usr/bin/env python
from PseudoNetCDF.plotutil.vertprofile import vertprofileplot, vertparser
from PseudoNetCDF.pncparse import pncparse

if __name__ == '__main__':
    ifiles, args = pncparse(has_ofile = True, parser = vertparser)
    vertprofileplot(ifiles, args)
