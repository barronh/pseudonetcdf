# PseudoNetCDF Installation #

## Known Platforms ##
  * Unix
  * Mac
  * Linux

## Requirements ##
**All requirements are met by using [Enthought Python Distribution](http://www.enthought.com/products/edudownload.php)**
Alternatively, you can install the following:
  1. [Python >=2.5](http://python.org)
  1. [NumPy >=1.2](http://numpy.scipy.org)
  1. [YAML](http://www.yaml.org)
  1. Optionally, one of the below netcdf readers
    * [pynetcdf](http://pypi.python.org/pypi/pynetcdf/0.7)
    * [pupynere](http://pypi.python.org/pypi/pupynere/)
    * [netCDF4](http://code.google.com/p/netcdf4-python)
    * [netCDF3](http://code.google.com/p/netcdf4-python)
    * [Scientific](http://dirac.cnrs-orleans.fr/plone/software/scientificpython/)
## Instructions ##
### General ###
# [Click here](https://code.google.com/p/pseudonetcdf/source/browse/) to browse source, at top click on "zip" beside download
  1. unzip PseudoNetCDF-src.zip -d PseudoNetCDF-src
  1. cd PseudoNetCDF-src/trunk
  1. **python** setup.py install

`*` **python** must be the executable with installed requirements

### Linux Without Permissions ###
  1. Navigate in terminal to the folder you want to install/work in.
  1. type `curl -kLO http://repo.continuum.io/miniconda/Miniconda-3.5.5-Linux-x86_64.sh`
    1. path may change in the future
  1. type `bash Miniconda-3.5.5-Linux-x86_64.sh`
    1. confirm prompts as necessary until you get to the install path (at the time of this writing: hit enter; hit space; type yes; hit enter;)
    1. Choose `./miniconda` as the install path
    1. Decide whether or not all future sessions should use this python
      1. Choose yes to enable this miniconda for all future sessions;
      1. Choose no to require miniconda to be enabled as needed;
  1. Enable miniconda (this step is required for at least this session): `export PATH=${PWD}/miniconda/bin:${PATH}`
  1. Install packages (type code in each step and confirm prompts)
    1. `conda install numpy`
    1. `conda install scipy`
    1. `conda install matplotlib`
    1. `conda install basemap`
    1. `conda install netcdf4`
    1. `conda install pip`
    1. `pip install --no-deps --upgrade --install-option="--prefix=${PWD}/miniconda/" git+https://code.google.com/p/pseudonetcdf/`
  1. You're ready to go!