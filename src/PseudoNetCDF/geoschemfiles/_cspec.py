import os
import re
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from PseudoNetCDF.geoschemfiles import bpch


def cspec(path, smv2='smv2.log'):
    if not os.path.isfile(smv2):
        if os.path.isdir(smv2):
            smv2 = os.path.join(smv2, 'smv2.log')
        else:
            smv2 = os.path.join(os.path.dirname(path), smv2)
    diaginfo, tracerinfo = get_info_for_cspec(get_cspec(smv2))
    return bpch(path, diaginfo=diaginfo, tracerinfo=tracerinfo)


def get_cspec(smvpath):
    text = open(smvpath, 'r').read()
    reo = re.compile(
        r'(?:NBR NAME\s+MW BKGAS\(VMRAT\))(?P<active>.*)' +
        r'(?:INACTIVE SPECIES FOR THIS RUN ARE:)' +
        r'(?P<inactive>.*)(?:THE DEAD SPECIES FOR THIS RUN ARE:)' +
        r'(?P<dead>.*?)(?:=====)', re.M | re.DOTALL)
    gv = reo.search(text).groupdict()
    active = [v.strip().split()[1] for v in gv['active'].strip().split('\n')]
    inactive = gv['inactive'].replace('\n', ' ').split()
    # dead = gv['dead'].replace('\n', ' ').split()
    all = active + inactive  # + dead
    return all


def get_info_for_cspec(all):
    out = ""
    for si, spc in enumerate(all):
        out += "%-8s %-8s cspec                 " % (spc, spc)
        out += "0.000E-00  1 %-8d 1.000E+00 molec/cm3\n" % (si + 1)

    initstr = ("#\n       0 IJ-CHK-$                                 " +
               "Tracer concentration                    \n")
    diag = StringIO(initstr)
    return diag, StringIO(out)


if __name__ == '__main__':
    specs = get_cspec('testdata/smv2.log')
    diagfile, trcfile = get_info_for_cspec(specs)
    f = bpch('/Users/barronh/Development/testrun/BC.CSPEC.20050101',
             tracerinfo=trcfile, diaginfo=diagfile)
