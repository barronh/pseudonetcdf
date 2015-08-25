import re
import unittest

#site-packages
import yaml
from numpy import array,zeros,newaxis
from PseudoNetCDF.sci_var import PseudoNetCDFFile as pncf, PseudoNetCDFVariable as pncv

net_num_re='(?<=\s)[-]?\d+[.]?\d+'
net_output_re=re.compile('^\w+(\s+'+net_num_re+'){24}$')
net_name_re=re.compile('[+]|source|termination|production|Source|Termination|Production|Oxidation|oxidation|Radical|radical')
net_time_dim_re=re.compile('^Ending Hour')
net_spc_name_re=re.compile('^\w+')
net_time_re=re.compile(net_num_re+'|daily|Daily')
net_num_re=re.compile(net_num_re)

sum_num_pattern='([-]?\d+[.]\d+)'
sum_num_sep='\s+'

def print_net_rxns(net_rxns,time='Daily'):
    nrxns=net_rxns.NET_RXN
    for nrxn in nrxns:
        print "%s: %s" % (nrxn,net_reaction(net_rxns,nrxn,time))

class sum_reader(pncf):
    @classmethod
    def isMine(cls, *args, **kwds):
        return False
    def __init__(self,sumfile):
        if type(sumfile)==str:
          self.sumfile=open(sumfile,'r')
        else:
          self.sumfile = sumfile
        if not 'Sum'.lower() in self.sumfile.readline().lower():
            raise IOError('Not a sum file')
        self.sumfile.seek(0, 0)
        time,var_names,values=self.parse()
    
        self.createDimension('TSTEP',len(time))
    
        for vn,vv in zip(var_names,values):
            tv=self.createVariable(vn,'f',('TSTEP',))
            tv[:] = vv

        tv=self.createVariable('TFLAG','S10',('TSTEP',))
        tv[:] = time
    
    def parse(self):
        time=[]
        var_names=[]
        values=[]
        lines=self.sumfile.readlines()
        for line in lines:
            line=line[:-1]
            if net_time_dim_re.match(line)!=None and time==[]:
                time=list(set(line.split(' ')[2:]))[1:]
                time.sort()
            sum_num_re=re.compile('(.*)\s+'+sum_num_sep.join([sum_num_pattern]*len(time)))

            if time!=[] and sum_num_re.match(line)==None and line!='':
                category=line
        
            if time!=[] and sum_num_re.match(line)!=None:
                grps=sum_num_re.match(line).groups()
                var_names.append(category.strip()+': '+grps[0].strip())
                values.append(map(float,grps[1:]))
  
        return time,var_names,values
      
def net_reaction(net_rxns,net_rxn,time='Daily'):
    spcs=net_rxns.SPECIES
    reactants=[]
    products=[]
    tmp="%f %s"
    for spc in spcs:
        v=net_rxns(NET_RXN=net_rxn,TIME=time,SPECIES=spc)
        if v==0:
            pass
        elif v>0:
            products.extend((v,spc))
        elif v<0:
            reactants.extend((-1*v,spc))
      
    return (" + ".join([tmp for i in reactants if type(i)==str])+" <-> "+" + ".join([tmp for i in products if type(i)==str])) % tuple(reactants + products)

class ctb_reader(pncf):
    @classmethod
    def isMine(cls, *args, **kwds):
        return False

    def __init__(self,file):
        if type(file)==str:
          file=open(file)
        lines=file.readlines()
        parameters,daily,peak=self.parse(lines)


        self.createDimension('PARAMETERS',len(parameters))
        self.createDimension('PEAK_DAILY',2)

        v=self.createVariable('CTB','f',('PARAMETERS','PEAK_DAILY'))
        v[:] = array([daily,peak]).swapaxes(0,1)
        v.units='ppbV'
    
    def parse(self,lines):
        lines=lines[16:]
        parameters,daily,hourly=[],[],[]
        for l in lines:
            try:
                p,d,h=l[:-1].split('\t')
                parameters.append(p)
                daily.append(float(d))
                hourly.append(float(h))
            except:
                break
        return parameters,daily,hourly


class net_reader(pncf):
    @classmethod
    def isMine(cls, *args, **kwds):
        return False

    def __init__(self,infile):
        self.spc=[]
        self.nrxns=[]
        self.time=[]

        netdict={}
        if type(infile)==str:
            infile=open(infile)
        lines=infile.readlines()
        self.parse(lines,False)
        self.createDimension('NET_RXN',len(self.nrxns))
        self.createDimension('SPECIES',len(self.spc))
        self.createDimension('TSTEP',len(self.time))
        self.createVariable('NET','f',('TSTEP','NET_RXN','SPECIES'))
        self.parse(lines,True)

    def parse(self,lines,fill=False):
        for l in lines:
          if net_time_dim_re.match(l)!=None and self.time==[]:
              self.time=net_time_re.findall(l)
          elif net_name_re.search(l)!=None:
              netrxn=l.strip()
              if netrxn not in self.nrxns:
                  self.nrxns.append(netrxn)
          elif net_output_re.match(l)!=None:
              spcn=net_spc_name_re.search(l).group()
              if spcn not in self.spc:
                  self.spc.append(spcn)
              if fill:
                  slc=self.variables['NET'][:,self.nrxns.index(netrxn),self.spc.index(spcn)]
                  slc[:]=map(float,net_num_re.findall(l))
                  
class mrgaloft(pncf):
    import re
    # create regular expressions to read the input file..
    sci_not = '-?\d+.\d+[E,e][+,-]\d+'
    time_re  = re.compile('Time =([0-9][0-9]0000)', re.IGNORECASE)
    irr_re   = re.compile('\{\s*(\d+)\}\s+('+sci_not+')', re.IGNORECASE)
    ipr_re   = re.compile('"\w+\s*"\s*('+sci_not+'\s+)+', re.IGNORECASE)
    split_on_blanks_re = re.compile('[ ]+')

    @classmethod
    def isMine(cls, *args, **kwds):
        return False

    def __init__(self,mrgfile):
        if isinstance(mrgfile,str):
            mrgfile=open(mrgfile,'rb')
        elif isinstance(mrgfile,file):
            pass
        else:
            raise TypeError, "mrgfile must be a path or a file"

        self.mrgfile=mrgfile
        self.__initfile()
        self.__get_spc_prc()
        self.__internalize()
        self.__initfile()
    
    def __get_spc_prc(self):
        #Read and discard time and irr
        self.read_time()
        self.read_irr()
    
        #Read ipr with species set to true
        #discard process information and internalize spc names
        self.spc,self.prc=self.read_ipr(True,True)[1:]
        self.Species="".join([i.ljust(16) for i in self.spc])
        self.Process="".join([i.ljust(16) for i in self.prc])

    def __initfile(self):
        #Move to start of file
        self.mrgfile.seek(0)
        #Internalize ipr and irr file names
        self.iprfilename=self.readline()[:-1]
        self.irrfilename=self.readline()[:-1]
    
    def __internalize(self):
        self.times=[]
        self.irr=[]
        self.ipr=[]
        #Internalize all arrays
        for t,r,p in self:
            self.times.append(t)
            self.irr.append(r)
            self.ipr.append(p)
        self.createDimension('TSTEP',len(self.times))
        self.createDimension('RXNS',len(self.irr[0]))
        self.createDimension('VAR',2)
        self.createDimension('DATE-TIME',2)
        self.createDimension('SPECIES',len(self.ipr[0]))
        self.createDimension('PROCESS',len(self.ipr[0][0]))
        v=self.createVariable('IRR','f',('TSTEP','RXNS'))
        v[:] = array(self.irr)
        v.units='ppb'
        v=self.createVariable('IPR','f',('TSTEP','SPECIES','PROCESS'))
        v[:] = array(self.ipr)
        v.units='ppb'
        v=self.createVariable('TFLAG','I',('TSTEP','VAR','DATE-TIME'))
        v[:,:,:]=array((zeros(len(self.times),'i'),self.times),'i').swapaxes(0,1)[:,newaxis,:]
        v[:,:,1]

    def readline(self):
        return self.mrgfile.readline()
  
    def read_time(self):
        f=self.mrgfile
        line = f.readline()
        if line[0].strip() == '|':
            raise EOFError('EOF File: found |')
    
        # next line is first Time =
        ir_time=self.time_re.match(line).groups()[0]
        if ir_time == None:
            raise ValueError, "ERROR:: in get_data read, did not find a time."
        return int(ir_time)
    
    def read_irr(self):
        f=self.mrgfile
        # initialize ir and ip storage...
        ir_rates = []  # account for 1-origin of rates
        
        # read in the !"Rxn No"    "Int Rate" header
        line = f.readline()
    
        line = f.readline()
    
        # read the irr values
        while line[0] != ';' :
            ir_value= self.irr_re.match(line).groups()[1]
            if ir_value == None:    # if this is a { n} n.nnnn line ...
                raise ValueError, "Expecting irr formatted line (i.e. { n} n.nnnn line ...)"
            else:
                ir_rates.append(float(ir_value))
            line = f.readline()
        return ir_rates
    
    def read_ipr(self,spc=False, prc=False):
        f=self.mrgfile
        # read in the ! Species      Initial conc. ... header
        line = f.readline()[15:]
        r=range(len(line)/17)
        ip_prc=[line[i*17:i*17+17].strip() for i in r]
    
        # 
        # ip_rates = [ [ip_values], [ip_values], [ip_values], ... ]
        #   order of sub-lists 'ip_values' is order of iProcess id Number
        #   for this hour.
        line = f.readline()
        ip_rates=[]
        ip_spc=[]
        while line[0] != ';' :
            ir_set=self.ipr_re.match(line)
            if ir_set != None:
                ip_spc.append(line[1:9].strip())
                ip_set = self.split_on_blanks_re.split(line[11:-1].strip())
                ip_values = [ float(v) for v in ip_set ]
                ip_rates.append(ip_values)
            line = f.readline()
        if spc and prc:
            return ip_rates,ip_spc,ip_prc
        elif spc:
            return ip_rates,ip_spc
        elif prc:
            return ip_rates,ip_prc
        else:
            return ip_rates

    def __iter__(self):
        self.__initfile()
        while True:
            try:
                yield self.read_time(),self.read_irr(),self.read_ipr()
            except EOFError:
                raise StopIteration 

class TestReaders(unittest.TestCase):
    def setUp(self):
        from PseudoNetCDF.testcase import net_balance_paths
        self.netfile=net_balance_paths['net_file']
        self.mrgfile=net_balance_paths['mrg_file']
        
    def testNet(self):
        netfile=net_reader(open(self.netfile))
        
    def testMrg(self):
        mrgfile=mrgaloft(open(self.mrgfile))
        irr=mrgfile.variables['IRR']
        ipr=mrgfile.variables['IPR']
        self.assertEqual(irr[0,60],array(2.34061E-02,'f')[...])
        self.assertEqual(irr[6,60],array(5.21899E-04,'f')[...])

        self.assertEqual(ipr[0,6,3],array(1.290572E+01,'f')[...])
        self.assertEqual(ipr[6,6,3],array(1.497632E-01,'f')[...])
        
    def runTest(self):
        pass
        
if __name__ == '__main__':
    unittest.main()
    
