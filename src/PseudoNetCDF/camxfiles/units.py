from PseudoNetCDF.camxfiles.aerosol_names import aerosol_names

_camx_units = {'EMISSIONS ': {True: 'g/time', False: 'mol/time'},
               'AVERAGE   ': {True: 'micrograms/m**3', False: 'ppm'},
               'BOUNDARY  ': {True: 'micrograms/m**3', False: 'ppm'},
               'INSTANT   ': {True: 'micrograms/m**3', False: 'micromoles/m**3'},
               'IPR': {True: 'micrograms/m**3', False: 'micromoles/m**3'},
               'AIRQUALITY': {True: 'micrograms/m**3', False: 'micromoles/m**3'}}

def get_uamiv_units(filename, key):
    from PseudoNetCDF.camxfiles.aerosol_names import aerosol_names
    return _camx_units[filename][key in aerosol_names]

