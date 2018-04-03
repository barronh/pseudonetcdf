from collections import defaultdict

__doc__ = r"""
.. _units
:mod:`units` -- Functions for converting units
========================================================

.. module:: units
   :platform: Unix, Windows
   :synopsis: Provides functions for converting units.
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""

__all__ = ['F2C', 'F2K', 'K2C', 'K2F', 'KCMAQ2F', 'M2km', 'MPS2kph',
           'convert', 'converter', 'converters_dict', 'km2m', 'm2ft',
           'm2km', 'm2miles', 'min2h', 'molespsCMAQ2molesph', 'mps2kmps',
           'mps2kph', 'mps2milesph', 'mps2milesps', 's2h', 's2min']


def s2min(a):
    return a / 60.


def min2h(a):
    return a / 60.


def s2h(a):
    return min2h(s2min(a))


def m2km(a):
    return a * .001


def M2km(a):
    return a * .001


def km2m(a):
    return a * 1000.


def m2ft(a):
    return a * 3.2808399


def m2miles(a):
    return a * 0.000621371192


def K2C(a):
    return a - 273.2


def K2F(a):
    return K2C(a) * 9 / 5. + 32


def KCMAQ2F(a):
    return K2C(a) * 9 / 5. + 32


def F2C(a):
    return (a - 32.) * 5 / 9.


def F2K(a):
    return F2C(a) + 273.2


def mps2kph(a):
    return a * 3.6


def MPS2kph(a):
    return a * 3.6


def mps2kmps(a):
    return a * .001


def mps2milesph(a):
    return a * 2.23693629


def mps2milesps(a):
    return a * 0.000621371192


def molespsCMAQ2molesph(a):
    return a * 60 * 60


def ppm2ppb(a):
    return a * 1000.


def ppb2ppt(a):
    return a * 1000.


def ppt2ppb(a):
    return a / 1000.


def ppb2ppm(a):
    return a / 1000.


def ppm2ppt(a):
    return ppm2ppb(ppb2ppt(a))


def ppt2ppm(a):
    return ppt2ppb(ppb2ppm(a))


class converters_dict(defaultdict):
    def __init__(self, dct):
        for k, v in dct.items():
            self[k] = v

    def __missing__(self, key):
        if type(key) == tuple and key[0] == key[1]:
            return lambda a: a


converter = converters_dict({
    ('s', 'min'): s2min,
    ('min', 'h'): min2h,
    ('s', 'h'): s2h,
    ('m/s', 'km/h'): mps2kph,
    ('M/S             ', 'km/h'): MPS2kph,
    ('m/s', 'miles/h'): mps2milesph,
    ('m/s', 'miles/s'): mps2milesps,
    ('m', 'km'): m2km,
    ('M               ', 'km'): M2km,
    ('km', 'm'): km2m,
    ('m', 'ft'): m2ft,
    ('m', 'miles'): m2miles,
    ('moles/s         ', 'moles/h'): molespsCMAQ2molesph,
    ('F', 'K'): F2K,
    ('F', 'C'): F2C,
    ('K', 'C'): K2C,
    ('K', 'F'): K2F,
    ('K               ', 'deg_F'): KCMAQ2F,
    ('K', 'deg_F'): K2F,
    ('ppm', 'ppb'): ppm2ppb,
    ('ppb', 'ppt'): ppb2ppt,
    ('ppt', 'ppb'): ppt2ppb,
    ('ppb', 'ppm'): ppb2ppm,
    ('ppt', 'ppm'): ppt2ppm,
    ('ppm', 'ppt'): ppm2ppt})


def convert(var, inunit, outunit):
    """
    Function that converts var from inunit to outunit.

    both units must be in converter dictionary and be compatible
    """
    return converter[(inunit, outunit)](var)
