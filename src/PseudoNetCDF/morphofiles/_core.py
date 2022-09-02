#!/usr/bin/env python -i
from __future__ import print_function

import re
from warnings import warn
from datetime import datetime

from PseudoNetCDF import PseudoNetCDFFile

#  20 lines reads hourly irr

_split1 = 'Environment Tables for'


def MorphoIRRt(irrpath):
    mrglines = open(irrpath).readlines()
    try:
        datelines = [_l for _l in mrglines if _split1 in _l]
        datestr = datelines[0].split(_split1)[-1].strip()
        jday = int(datetime.strptime(datestr, '%d-%b-%y').strftime('%Y%j'))
    except Exception:
        warn('Could not find/parse date; using 1900001')
        jday = 1900001

    mrglines = [line for line in mrglines if line[:2]
                not in ('//', '**') and line not in ('', '\n')]
    name_line = mrglines.pop(0)
    irrlabel = re.compile(r'rt\[\S+\]')
    rxn_names = ['IRR_%s' % name.replace('rt[', '').replace(']', '')
                 for name in irrlabel.findall(name_line)]
    name_line = ['N', 'T'] + rxn_names
    unit_line = mrglines.pop(0).split()
    unit_line = [unit for unit in unit_line]
    assert (all([eval(v) == 0. for v in mrglines.pop(0).split()][2:]))
    mrgfile = PseudoNetCDFFile()
    mrgfile.createDimension('TSTEP', len(mrglines))
    mrgfile.createDimension('DATE-TIME', 2)
    mrgfile.createDimension('VAR', 1)
    unit_dict = dict([(k, v) for k, v in zip(name_line, unit_line)])
    tflag = mrgfile.createVariable('TFLAG', 'f', ('TSTEP', 'VAR', 'DATE-TIME'))
    tflag.units = '<JDAY, MIN>'
    tflag.long_name = tflag.var_desc = 'TFLAG'
    tflag[:, :, 0] = jday
    for name in name_line:
        var = mrgfile.createVariable(name, 'f', ('TSTEP',))
        var.units = unit_dict[name]
        var.long_name = var.var_desc = name
    for ti, line in enumerate(mrglines):
        for var_name, value in zip(name_line, line.split()):
            var = mrgfile.variables[var_name]
            var[ti] = float(value)

    for name in name_line:
        if name in ('T', 'N'):
            continue
        var = mrgfile.variables[name]
        var[1:] = (var[1:] - var[:-1]) * 1000.
    tflag[:, :, 1] = mrgfile.variables['T'][:, None]
    return mrgfile


def MorphoConc(concpath):
    conclines = open(concpath).readlines()
    conclines = [line for line in conclines if line[:2]
                 not in ('//', '**') and line not in ('', '\n')]
    name_line = conclines.pop(0)
    conclabel = re.compile(r'n\[\S+\]')
    conc_names = [name.replace('n[', '').replace(']', '')
                  for name in conclabel.findall(name_line)]
    name_line = ['N', 'T'] + conc_names
    unit_line = conclines.pop(0).split()
    unit_line = [unit for unit in unit_line]
    concfile = PseudoNetCDFFile()
    concfile.createDimension('TSTEP', len(conclines))
    concfile.createDimension('DATE-TIME', 2)
    concfile.createDimension('VAR', 1)
    unit_dict = dict([(k, v) for k, v in zip(name_line, unit_line)])
    for name in name_line:
        var = concfile.createVariable(name, 'f', ('TSTEP',))
        var.units = unit_dict[name]
        var.long_name = var.var_desc = name
    for ti, line in enumerate(conclines):
        for var_name, value in zip(name_line, line.split()):
            var = concfile.variables[var_name]
            var[ti] = float(value)

    return concfile


def MorphoPERMM(concpath, irrtpath):
    mrgf = MorphoIRRt(irrtpath)
    concf = MorphoConc(concpath)
    for key, var in concf.variables.items():
        initvar = mrgf.createVariable('INIT_%s' % key, 'f', ('TSTEP',))
        initvar.long_name = var.long_name
        initvar.var_desc = var.var_desc
        initvar.units = var.units
        initvar[:] = var[:-1]
        finalvar = mrgf.createVariable('FINAL_%s' % key, 'f', ('TSTEP',))
        finalvar.long_name = var.long_name
        finalvar.var_desc = var.var_desc
        finalvar.units = var.units
        finalvar[:] = var[1:]
    mrgf.Species = ' '.join([k.ljust(16) for k in concf.variables.keys()])
    mrgf.Processes = 'INIT'.ljust(17) + 'FINAL'.ljust(16)
    return mrgf


# if __name__ == '__main__':
#     from permm.guis.Simplewx import StartWx
#     import os
#     from permm.analyses.network.core import MakeCarbonTrace
#     from permm import Mechanism
#     concpath = sys.argv[1]
#     mrgpath = sys.argv[2]
#     irrtpath = sys.argv[3]
#     print(sys.argv[1:])
#     if os.path.exists('mech.yaml'):
#         cb05 = Mechanism('mech.yaml')
#     else:
#         cb05 = MorphoMrg(mrgpath, 'mech.yaml')
#     mrg = MorphoPERMM(concpath, irrtpath)
#     cb05.set_mrg(mrg)
#     cb05.globalize(globals())
#     n = MakeCarbonTrace(cb05, makes_larger=[PAR])
#     import networkx as nx
#     es = nx.dfs_edges(n)
#     #from permm.analyses.history import matrix
#     #history = matrix(cb05, [C2O3], [HC+Radical-OH-HO2-O1D], [])
