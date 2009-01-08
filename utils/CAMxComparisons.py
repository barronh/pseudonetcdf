"""
.. _CAMxComparisons
:mod:`CAMxComparisons` -- CAMx Comparisons utilities
====================================================

.. module:: CAMxComparisons
   :platform: Unix, Windows
   :synopsis: Provides CAMx Comparison utilities used in :ref:`pyPA`
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""
__all__ = ['grid_diff', 'grid_diff_sum', 'point_compare', 'emiss_compare']

HeadURL="$HeadURL$"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy$"
__version__ = RevisionNum

from pyPA.utils.CAMxFiles import gridded_emissions,point_source
from pyPA.utils.sci_var import PseudoNetCDFFile as pncf

def grid_diff(outfile,infile1,infile2):
    """Grid compare compares gridded emission files
    using file1-file 2 and returns either an hourly
    or daily sum.  Species moles are converted to
    mass using a species weight factor
    file1 = path string or open file like object
    file2 = path string or open file like object
    spc_wt = array dim(species)
    hourly = boolean true=hourly false=daily
    """
    g1=gridded_emissions(infile1)
    g2=gridded_emissions(infile2)
    outfile=Pseudo2NetCDF().convert(g2,outfile)
    for k,v in g1.variables.iteritems():
        outfile.variables[k]-=g2.variables[k]
    return outfile

def grid_diff_sum(outfile,infile1,infile2,hourly=True):
    tmppncf=PseudoNetCDFFile()
    tmppncf=grid_diff(tmppncf,infile1,infile2)
    Pseudo2NetCDF().addGlobalProperties(tmpncf,outfile)
    Pseudo2NetCDF().addDimensions(tmpncf,outfile)
    for k,v in tempncf.variables.iteritems():
        tv=outfile.createVariable(k,'f',{True: ('TSTEP',),False: ()}[hourly])
        if hourly:
            tv.assignValue(v.sum(-1).sum(-1).sum(-1))
        else:
            tv.assignValue(v.sum())

def point_compare(outfile,infile1,infile2,spc_wt=1,hourly=True):
    """Point compare compares point emission files
    using file1-file 2 and returns either an hourly
    or daily sum.  Species moles are converted to
    mass using a species weight factor
    file1 = path string or open file like object
    file2 = path string or open file like object
    spc_wt = array dim(species)
    hourly = boolean true=hourly false=daily
    """
    p1=point_source(infile1)
    p2=point_source(infile2)
    pa1=p1.getArray()
    pa2=p2.getArray()
    mol_diff=pa1-pa2
    mol_diff=mol_diff.sum(axis=-1)
    if not hourly:
        mol_diff=mol_diff.sum(0)
    return mol_diff*spc_wt

def emiss_compare(gfile1,pfile1,gfile2,pfile2,spc_wt=1,hourly=True):
    """Emiss compare compares emission files
    using file1-file 2 and returns either an hourly
    or daily sum.  Species moles are converted to
    mass using a species weight factor
    gfile1 = path string or open file like object
    gfile2 = path string or open file like object
    pfile1 = path string or open file like object
    pfile2 = path string or open file like object
    spc_wt = array dim(species)
    hourly = boolean true=hourly false=daily
    """
    p=point_compare(pfile1,pfile2,spc_wt,hourly)
    g=grid_compare(gfile1,gfile2,spc_wt,hourly)
    return p+g
