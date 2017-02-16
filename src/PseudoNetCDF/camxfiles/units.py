__all__ = ['get_chemparam_names', 'get_uamiv_units']

_camx_units = {'EMISSIONS ': {True: 'g/time', False: 'mol/time'},
               'AVERAGE   ': {True: 'micrograms/m**3', False: 'ppm'},
               'BOUNDARY  ': {True: 'micrograms/m**3', False: 'ppm'},
               'INSTANT   ': {True: 'micrograms/m**3', False: 'micromoles/m**3'},
               'IPR': {True: 'micrograms/m**3', False: 'micromoles/m**3'},
               'AIRQUALITY': {True: 'micrograms/m**3', False: 'ppm'},
               'DEPOSITION': {True: 'g/ha', False: 'mol/ha'},}

def get_chemparam_names(chemparampath):
    inlines = open(chemparampath, 'rU').readlines()
    startaero = None
    for li, l in enumerate(inlines):
        if l[:14] == "     Gas Spec ":
            startgas = li + 1
        if l[:14] == "     Aero Spec":
            endgas = li
            startaero = li + 1
        if l[:16] == "Reaction Records":
            if startaero is None:
                endgas = li
            endaero = li
            break
    
    gaskeys = tuple([l[5:15].strip() for l in inlines[startgas:endgas]])
    if startaero is None:
        return dict(gas = gaskeys, aer = ())
    aerkeys = tuple([l[5:15].strip() for l in inlines[startaero:endaero]])
    return dict(gas = gaskeys, aerosol = aerkeys)

def get_uamiv_units(filename, key, aerosol_names = None):
    if aerosol_names is None:
        from PseudoNetCDF.camxfiles.aerosol_names import aerosol_names
    if key[-3:] in ('_DV',):
        return 'm/s'
    elif key[-3:] in ('_DD', '_WD'):
        return _camx_units['DEPOSITION'][key[:-3] in aerosol_names]
    return _camx_units[filename][key in aerosol_names]

