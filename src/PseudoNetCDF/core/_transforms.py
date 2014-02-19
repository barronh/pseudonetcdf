from _files import PseudoNetCDFFile
from _variables import PseudoNetCDFVariable
from PseudoNetCDF.units import convert


def PseudoNetCDFVariableConvertUnit(var,outunit):
    """
    Convert the unit of var and update the 
    associated IOAPI metadata
    """
    do = PseudoNetCDFFile()
    shape=var.shape
    for i,d in enumerate(var.dimensions):
        do.createDimension(d, shape[i])
    outvar=PseudoNetCDFVariable(do,var.long_name.strip(),var.typecode(),var.dimensions,values=convert(var,var.units,outunit))
    for k in var.ncattrs():
        v = getattr(var, k)
        setattr(outvar,k,v)
    outvar.units=outunit
    return outvar
    
