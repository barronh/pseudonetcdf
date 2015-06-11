try:
    from setuptools import setup
except:
    from distutils.core import setup
import os
import sys
from warnings import warn
netcdf_pkgs = [('netCDF4', 'Dataset', 'Variable'), \
               ('netCDF3', 'Dataset', 'Variable'), \
               ('pupynere', 'NetCDFFile', 'NetCDFVariable')]
for pkg, reader, variable in netcdf_pkgs:
    try:
        NetCDFFile = getattr(__import__(pkg, fromlist = [reader]),reader)
        define_function = "from %s import %s as NetCDFFile" % (pkg, reader)
        define_variable = "from %s import %s as NetCDFVariable" % (pkg, variable)
        netcdfpkg = [pkg]
        break
    except ImportError, e:
        warn(e.message)
else:
    warn("Did not find a 'true' NetCDFFile reader; 'true' NetCDF functionality will be disabled")
    netcdfpkg = []
    define_function = """
class NetCDFFile(object):
    def __init__(self, *args, **kwds):
        raise ImportError('System has no valid netCDF reader; install netcdf4-python or pupynere')
"""

print >> file(os.path.join('src', 'PseudoNetCDF', 'netcdf.py'),'wb'), """
__all__ = ['NetCDFFile']
__doc__ = \"\"\"
.. _netcdf
:mod:`netcdf` -- netcdf import point
====================================

.. module:: netcdf
   :platform: Unix, Windows
   :synopsis: Povides a single import point for a package.  If
              a user has one of many netcdf interfaces, this module
              selects it and provides it.
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
\"\"\"
%s
%s
""" % (define_function, define_variable)

def find_packages():
    import os
    packages = []
    walker = os.walk('src')
    prefix = os.path.join(os.path.curdir,'src')
    for thisdir, itsdirs, itsfiles in walker:
        if '__init__.py' in itsfiles:
            packages.append(thisdir[len(prefix)-1:])
    
    return packages
            
def find_data():
    import os
    import re
    data_pattern = re.compile(r'.*(.|_)(yaml|nc|net|irr|phy|ptb|sum|voc|txt|xls|graffle)$')
    data = []
    prefix = os.path.join(os.path.curdir,'src', 'PseudoNetCDF')
    walker = os.walk('src')
    for thisdir, itsdirs, itsfiles in walker:
        if thisdir != os.path.join('src','PseudoNetCDF.egg-info'):
            data.extend([os.path.join(thisdir[len(prefix)-1:],f) for f in itsfiles if data_pattern.match(f) is not None])
    
    return data

packages = find_packages()
data = find_data()


setup(name = 'PseudoNetCDF',
      version = '1.0',
      author = 'Barron Henderson',
      author_email = 'barronh@gmail.com',
      maintainer = 'Barron Henderson',
      maintainer_email = 'barronh@gmail.com',
      packages = packages,
      package_dir = {'': 'src'},
      package_data = {'PseudoNetCDF': data},
      scripts = ['scripts/aqsraw4pnceval.py', 'scripts/aqsrest4pnceval.py', 'scripts/asos4pnceval.py', 'scripts/pncboundaries.py', 'scripts/pncdump', 'scripts/pnceval', 'scripts/pncgen', 'scripts/pncload', 'scripts/pncmap.py', 'scripts/pncvertprofile.py', 'scripts/pncview'],
      requires = netcdfpkg + ['numpy (>=1.2)', 'yaml'],
      url = 'https://code.google.com/p/pseudonetcdf'
      )
