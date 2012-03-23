__all__=['osat']
__doc__ = """
.. _Write
:mod:`Write` -- CAMx uamiv variable transformations
==================================================

.. module:: Write
   :platform: Unix, Windows
   :synopsis: Provides :ref:`PseudoNetCDF` variable transformations
              for CAMx uamiv files.  See PseudoNetCDF.sci_var.PseudoNetCDFFile 
              for interface details
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""


from numpy import array,zeros
from warnings import warn
import unittest
from PseudoNetCDF.sci_var import PseudoNetCDFFile, PseudoNetCDFVariables, PseudoIOAPIVariable
from Memmap import uamiv

class osat(PseudoNetCDFFile):
    """
    OSAT provides a way of slicing and summing tagged species.
    
    sources - dictionary of name and id for CAMx source codes 
              (i.e. {'SOURCES1-5': '001','002','003','004','005'})
    regions - dictionary of name and id for CAMx region codes  
              (i.e. {'REGIONS1-5': '001','002','003','004','005'})
    
    Example:
        # Source and regions definitions
        sources={'SOURCES1-5': ('001','002','003','004','005')}
        regions={'REGIONS1-5': ('001','002','003','004','005')}
        # File initiation
        osatfile=osat('path',sources,regions)
        # Variable (v) would return the sum of all O3V for all five 
        # sources and all five regions
        v=osatfile.variables['O3V_SOURCES1-5_REGIONS1-5'] 
        
    """
    __delim='_'
    def __init__(self,rffile,sources={},regions={}):
        self.__child=uamiv(rffile)
        
        # Use Use all species based keys
        self.__sourcesbyNm=dict([(k[3:6],k[3:6]) for k in self.__child.variables.keys() if k[:3] in ('NOX','VOC','O3V','O3N')])
        self.__regionsbyNm=dict([(k[6:],k[6:]) for k in self.__child.variables.keys() if k[:3] in ('NOX','VOC','O3V','O3N')])
        
        # Create global keys
        self.__sourcesbyNm['']=tuple(self.__sourcesbyNm.keys())
        self.__regionsbyNm['']=tuple(self.__regionsbyNm.keys())
        
        for k,v in sources.iteritems():
            if type(v)==str:
                sources[k]=(v,)
        for k,v in regions.iteritems():
            if type(v)==str:
                regions[k]=(v,)
        
        # Update with user supplied keys
        self.__sourcesbyNm.update(sources)
        self.__regionsbyNm.update(regions)

        self.dimensions=self.__child.dimensions.copy()

        spc_keys=['NOX','VOC','O3N','O3V']
        time_keys=['I','D']
        time_dim_keys=[(k[:1],k[1:6],k[6:]) for k in time_keys]
        allkeys=[i for i in self.__child.variables.keys()]
        for skey in spc_keys:
            for src in sources.keys()+['']:
                for reg in regions.keys()+['']:        
                    allkeys.append(self.__delim.join([skey,src,reg]))
        allkeys=[i for i in set(allkeys)]
        self.variables=PseudoNetCDFVariables(self.__variables,allkeys)
    def __indices(self,keys):
        return [self.__child.__var_names__.index(k) for k in keys]

    def __var_id_names(self,var_name):
        components=var_name.split(self.__delim)
        if len(components)==1:
            keys=components
        else:
            spc,source,region=components
            keys=[]
            sources=self.__sourcesbyNm[source]
            regions=self.__regionsbyNm[region]
            for source in sources:
                for region in regions:
                    key=spc+source+region
                    if key in self.__child.variables.keys():
                        keys.append(key)
        return keys
        
    def __variables(self,key):
        var_id_names=self.__var_id_names(key)
        #var_nm_names=self.__var_nm_names(key)
        
        var_indices=self.__indices(var_id_names)
        cnt=False
        if len(var_id_names)>1:
            outvals=zeros(map(len, (self.dimensions['TSTEP'],self.dimensions['LAY'],self.dimensions['ROW'],self.dimensions['COL'])),dtype='>f')
            for k in var_id_names:
                outvals[...]+=self.__child.variables[k]
        else:
            outvals=self.__child.variables[var_id_names[0]]
            
        dimensions=('TSTEP','VAR','LAY','ROW','COL')
        v=PseudoIOAPIVariable(self,key,'f',dimensions,values=outvals, units = 'ppm')
        v.VAR_NAMES=''.join([nm.ljust(16) for nm in var_id_names])
        #v.VAR_NAME_DESCS=''.join([nm.ljust(16) for nm in var_nm_names])
        return v

class TestMemmap(unittest.TestCase):
    def runTest(self):
        pass
    def setUp(self):
        pass
    def testOSAT(self):
        from pyPA.testcase import CAMxOSATInst
        warn("OSAT test is not implemented: review output for reasonabilty")
        o=osat(CAMxOSATInst)
        k='NOX__'
        
        v=o.variables[k]
        print k,v.VAR_NAMES
        print v.sum(),v.min(),v.mean(),v.max()
        keys=[k for k in o._osat__child.variables.keys() if k[:3]=='NOX']
        for k in keys:
            v=o.variables[k]
            print k,v.sum(),v.min(),v.mean(),v.max()
