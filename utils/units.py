HeadURL="$HeadURL$"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy$"
__version__ = RevisionNum

try:
    from collections import defaultdict
except:
    from pyPA.utils.util import defaultdict

__all__ = ['F2C', 'F2K', 'K2C', 'K2F', 'KCMAQ2F', 'M2km', 'MPS2kph', 'convert', 'converter', 'converters_dict', 'km2m', 'm2ft', 'm2km', 'm2miles', 'min2h', 'molespsCMAQ2molesph', 'mps2kmps', 'mps2kph', 'mps2milesph', 'mps2milesps', 's2h', 's2min']

s2min=lambda a: a/60.
min2h=lambda a: a/60.
s2h=lambda a: min2h(s2min(a))
m2km=lambda a: a*.001
M2km=lambda a: a*.001
km2m=lambda a: a*1000.
m2ft=lambda a: a*3.2808399
m2miles=lambda a: a*0.000621371192
K2C=lambda a: a-273.2
K2F=lambda a: K2C(a)*9/5.+32
KCMAQ2F=lambda a: K2C(a)*9/5.+32
F2C=lambda a: (a-32.)*5/9.
F2K=lambda a: F2C(a)+273.2
mps2kph=lambda a: a*3.6
MPS2kph=lambda a: a*3.6
mps2kmps=lambda a: a*.001
mps2milesph=lambda a: a*2.23693629
mps2milesps=lambda a: a*0.000621371192
molespsCMAQ2molesph=lambda a: a*60*60

class converters_dict(defaultdict):
    def __init__(self,dct):
        for k,v in dct.iteritems():
            self[k]=v
    def __missing__(self,key):
        if type(key)==tuple and key[0]==key[1]:
            return lambda a: a
            
converter=converters_dict({('s','min'): s2min, \
                           ('min','h'): min2h, \
                           ('s','h'): s2h, \
                           ('m/s','km/h'): mps2kph, \
                           ('M/S             ','km/h'): MPS2kph, \
                           ('m/s','miles/h'): mps2milesph, \
                           ('m/s','miles/s'): mps2milesps, \
                           ('m','km'): m2km, \
                           ('M               ','km'): M2km, \
                           ('km','m'): km2m, \
                           ('m','ft'): m2ft, \
                           ('m','miles'): m2miles, \
                           ('moles/s         ','moles/h'): molespsCMAQ2molesph, \
                           ('F','K'): F2K, \
                           ('F','C'): F2C, \
                           ('K','C'): K2C, \
                           ('K','F'): K2F, \
                           ('K               ','deg_F'): KCMAQ2F, \
                           ('K','deg_F'): K2F} \
                           )

def convert(var, inunit, outunit):
    return converter[(inunit,outunit)](var)