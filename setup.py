from __future__ import print_function
try:
    from setuptools import setup
except Exception:
    from distutils.core import setup
import os
import subprocess

join = os.path.join


def find_packages():
    import os
    packages = []
    walker = os.walk('src')
    prefix = join(os.path.curdir, 'src')
    for thisdir, itsdirs, itsfiles in walker:
        if '__init__.py' in itsfiles:
            packages.append(thisdir[len(prefix) - 1:])

    return packages


def find_data():
    import os
    import re
    suffixes = ['yaml', 'nc', 'net', 'irr', 'phy', 'ptb',
                'sum', 'voc', 'txt', 'xls', 'graffle']
    data_pattern = re.compile(
        r'.*(.|_)(' + '|'.join(suffixes) + ')$')
    data = []
    prefix = join(os.path.curdir, 'src', 'PseudoNetCDF')
    walker = os.walk('src')
    for thisdir, itsdirs, itsfiles in walker:
        if thisdir != os.path.join('src',
                                   'PseudoNetCDF.egg-info'):
            data.extend([join(thisdir[len(prefix) - 1:], f)
                         for f in itsfiles
                         if data_pattern.match(f) is not None])

    return data


packages = find_packages()
data = find_data()
MAJOR = 3
MINOR = 1
MICRO = 0
ISRELEASED = False
VERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)


# Return the git revision as a string
def git_version():
    # Copied from scipy
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        Popen = subprocess.Popen
        out = Popen(cmd, stdout=subprocess.PIPE, env=env).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
        GIT_REVISION = out.strip().decode('ascii')
    except OSError:
        GIT_REVISION = "Unknown"

    return GIT_REVISION


def get_version_info():
    # Modified from scipy
    # Adding the git rev number needs to be done inside
    # write_version_py(), otherwise the import of
    # PseudoNetCDF.version messes up the build under Python 3.
    FULLVERSION = VERSION
    if os.path.exists('.git'):
        GIT_REVISION = git_version()
    elif os.path.exists('src/PseudoNetCDF/version.py'):
        # must be a source distribution, use existing version file
        # load as separate module to not load PseudoNetCDF/__init__.py
        import imp
        version = imp.load_source('PseudoNetCDF.version',
                                  'src/PseudoNetCDF/version.py')
        GIT_REVISION = version.git_revision
    else:
        GIT_REVISION = "Unknown"

    if not ISRELEASED:
        FULLVERSION += '.dev0+' + GIT_REVISION[:7]

    return FULLVERSION, GIT_REVISION


def write_version_py(filename='src/PseudoNetCDF/version.py'):
    # Modified from scipy
    cnt = """
# THIS FILE IS GENERATED FROM PseudoNetCDF SETUP.PY
short_version = '%(version)s'
version = '%(version)s'
full_version = '%(full_version)s'
git_revision = '%(git_revision)s'
release = %(isrelease)s
if not release:
    version = full_version
"""
    FULLVERSION, GIT_REVISION = get_version_info()

    a = open(filename, 'w')
    try:
        a.write(cnt % {'version': VERSION,
                       'full_version': FULLVERSION,
                       'git_revision': GIT_REVISION,
                       'isrelease': str(ISRELEASED)})
    finally:
        a.close()


write_version_py(filename='src/PseudoNetCDF/version.py')

long_desc = """NetCDF, NCO, and CDO are fantastic softwares that I use,
emulate, and admire. The primary drawbacks are restrictions on which
scientific data sources they will and won't work with, and what types
of operations they will and won't do. PseudoNetCDF was originally just
a NetCDF-like interface for many data formats, but has grown to include
many functionalities from NCO and CDO. This is a platform independent,
easy to install software to make scientific data manipulation easy."""

script_list = ['scripts/pncmadis2pnceval.py', 'scripts/pncaqsraw4pnceval.py',
               'scripts/pncaqsrest4pnceval.py', 'scripts/pncasos4pnceval.py',
               'scripts/pnc1d.py', 'scripts/pnc2d.py',
               'scripts/pncboundaries.py', 'scripts/pncdiurnal.py',
               'scripts/pncdump', 'scripts/pncdump.py', 'scripts/pnceval',
               'scripts/pnceval.py', 'scripts/pncgen', 'scripts/pncgen.py',
               'scripts/pncglobal2cmaq.py', 'scripts/pncload',
               'scripts/pncmap.py', 'scripts/pncqq.py',
               'scripts/pncscatter.py', 'scripts/pncts.py',
               'scripts/pncvertprofile.py', 'scripts/pncview',
               'scripts/pncview.py', 'scripts/pncwindrose.py']

requires_list = ['numpy>=1.2', 'netCDF4', 'scipy', 'matplotlib', 'pyyaml']
extra_requires_dict = {'textfiles': ['pandas'],
                       'projections': ['pyproj'],
                       'mapping': ['basemap'],
                       }

setup(name='PseudoNetCDF',
      version=VERSION,
      author='Barron Henderson',
      author_email='barronh@gmail.com',
      maintainer='Barron Henderson',
      maintainer_email='barronh@gmail.com',
      description=('Like NetCDF and NCO, but works with NetCDF and other ' +
                   'scientific formats.'),
      long_description=long_desc,
      packages=packages,
      package_dir={'': 'src'},
      package_data={'PseudoNetCDF': data},
      scripts=script_list,
      install_requires=requires_list,
      extras_require=extra_requires_dict,
      url='http://github.com/barronh/pseudonetcdf/',
      classifiers=['Programming Language :: Python :: 2.7',
                   'Programming Language :: Python :: 3',
                   'Operating System :: MacOS',
                   'Operating System :: Microsoft :: Windows',
                   'Operating System :: POSIX',
                   'Topic :: Scientific/Engineering',
                   'Topic :: Scientific/Engineering :: Atmospheric Science',
                   ]
      )
