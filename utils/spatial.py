HeadURL="$HeadURL$"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy$"
__version__ = RevisionNum

from math import cos,log,radians,degrees,pi,tan,sin,sqrt,atan
from warnings import warn

def pnpoly((x,y),poly):
    """
    pnpoly tests wether a point exists in a polygon or line
    The test is boundary inclusive.  In the case of a polygon, the
    last point should be the same as the first.
    """
    npol = len(poly)
    xp = [poly[i][0] for i in range(npol)]
    yp = [poly[i][1] for i in range(npol)]

#  in the special case of a line
    line = poly[0] != poly[-1]

    inpoly = False
    ptcount = range(npol)
    for i in range(npol):
        j = ptcount[i-1]
    
        if (y>=min(yp[i],yp[j]) and y<=max(yp[i],yp[j])):
            if yp[i] != yp[j]:
                xinters = (xp[j] - xp[i]) * (y - yp[i]) / (yp[j] - yp[i]) + xp[i]
                if x == xinters:
                    if not line:
                        inpoly = True
                        break
                    elif x>=min(xp[i],xp[j]) and x<=max(xp[i],xp[j]):
                        inpoly = True
                        break
                elif not line and x < xinters:
                    inpoly = not inpoly
            elif x>=min(xp[i],xp[j]) and x<=max(xp[i],xp[j]):
                inpoly = True
                break
                
    return inpoly

class poly(object):
    """Poly object is intended to provide pnpoly functionality
    (is point in polygon) without regenerating the polygon object
    each time.  It has been expanded for more functionality.
    """
    def __init__(self, poly):
        """Internalize the polygon object and save usefu
        attributes for future use.
        """
        
        self.poly = poly
        self.npol = len(poly)
        self.xp = [poly[i][0] for i in range(self.npol)]
        self.yp = [poly[i][1] for i in range(self.npol)]
        #  in the special case of a line
        self.line = poly[0] != poly[-1]
        self.ptcount = range(self.npol)

        self.points = self.cells()
        
    def cells(self):
        """Return a dictionary of all cells with
        -1 values.  Not sure if we are using the
        dictionary functionality anymore.
        """
        
        pts={}
        envelope = cartesian(range(min(self.xp),max(self.xp)+1),range(min(self.yp),max(self.yp)+1))
        while True:
            try:
                pt=envelope.next()
                if self.ptcheck(*pt):
                    pts[pt]=-1
            except:
                return pts

    def north(self):
        """Return list of cells on the top or north boundary
        """
        cells=self.cells()
        return [c for c in cells if not add_tuple(c,(0,1)) in cells]

    def south(self):
        """Return list of cells on the bottom or south boundary
        """

        cells=self.cells()
        return [c for c in cells if not add_tuple(c,(0,-1)) in cells]

    def east(self):
        """Return list of cells on the right or east boundary
        """

        cells=self.cells()
        return [c for c in cells if not add_tuple(c,(1,0)) in cells]

    def west(self):
        """Return list of cells on the left or west boundary
        """
        cells=self.cells()
        return [c for c in cells if not add_tuple(c,(-1,0)) in cells]

    def ptcheck(self,x,y):
        """Is this x,y inside the polygon property
        """
        inpoly = False
        for i in range(self.npol):
            j = self.ptcount[i-1]

            if (y>=min(self.yp[i],self.yp[j]) and y<=max(self.yp[i],self.yp[j])):
                if self.yp[i] != self.yp[j]:
                    xinters = (self.xp[j] - self.xp[i]) * (y - self.yp[i]) / (self.yp[j] - self.yp[i]) + self.xp[i]
                    if x == xinters:
                        if not self.line:
                            inpoly = True
                            break
                        elif x>=min(self.xp[i],self.xp[j]) and x<=max(self.xp[i],self.xp[j]):
                            inpoly = True
                            break
                    elif not self.line and x < xinters:
                        inpoly = not inpoly
                elif x>=min(self.xp[i],self.xp[j]) and x<=max(self.xp[i],self.xp[j]):
                    inpoly = True
                    break
                
        return inpoly

class gridconverter(object):
    def __init__(self,small,big):
        warn("gridconverter is deprecated favoring the DomainCoords class",DeprecationWarning)
        self.small=small
        self.big=big
        
    def ij(self,i,j):
        return self.big.ij(*self.small.xy(i,j))

class grid(object):
    def __init__(self,xorg,ncols,yorg,nrows,size,nlays,tops=None):
        warn("grid is deprecated favoring the DomainCoords class",DeprecationWarning)
        self.xorg=xorg
        self.yorg=yorg
        self.size=size
        self.ncols=ncols
        self.nrows=nrows
        self.nlays=nlays
        
    def x(self,i):
        if i<1 or i > self.ncols:
            raise KeyError, "Grid has only 1 thru %d cols" % self.ncols
        return (self.size*(i-1)+self.xorg)
        
    def y(self,j):
        if j<1 or j > self.nrows:
            raise KeyError, "Grid has only 1 thru %d rows" % self.nrows
        return (self.size*(j-1)+self.yorg)
        
    def xy(self,i,j):
        return self.x(i),self.y(j)
        
    def i(self,x):
        if x < self.xorg or x > self.x(self.ncols):
            return None
        else:
            return int((x-self.xorg)/self.size)+1

    def j(self,y):
        if y < self.yorg or y > self.y(self.nrows):
            return None
        else:
            return int((y-self.yorg)/self.size)+1
        
    def ij(self,x,y):
        return self.i(x), self.j(y)
        
class lccConv(object):
    """
    http://mathworld.wolfram.com/LambertConformalConicProjection.html
    """
    __slots__=['lon_ref','lat_ref','lat_std1','lat_std2','xy','latlon','F','N','rho','rho0','convin','EarthRadius','RHO0','sgnN','convout']
    def __init__(self,lon_ref=-100,lat_ref=40,lat_std1=30,lat_std2=60,lon=None,lat=None,x=None,y=None,EarthRadius=6370.,radiansFlg=False):
        warn("Deprecated favoring DomainCoords",DeprecationWarning)
        if radiansFlg:
            self.convin=lambda x: x
            self.convout=lambda x: x
        else:
            self.convin=radians
            self.convout=degrees
        
        self.lon_ref=self.convin(lon_ref)
        self.lat_ref=self.convin(lat_ref)
        self.lat_std1=self.convin(lat_std1)
        self.lat_std2=self.convin(lat_std2)
        self.EarthRadius=float(EarthRadius)
        
        self.N=self.n()
        self.sgnN=(-1,1)[self.n > 0]

        self.F=self.f()
        self.RHO0=self.rho0()
        
        if not (lon==None or lat==None):
            return self.xy(lat,lon)
        if not (x==None or y==None):
            return self.latlon(x,y)

    def xy(self,lat,lon):
        lat=self.convin(float(lat))
        lon=self.convin(float(lon))
        x=self.rho(lat)*sin(self.N*(lon-self.lon_ref))*self.EarthRadius
        y=(self.RHO0-self.rho(lat)*cos(self.N*(lon-self.lon_ref)))*self.EarthRadius
        return x,y
    
    def latlon(self,x,y,units=1000.):
        x=float(x)/self.EarthRadius
        y=float(y)/self.EarthRadius
        lat=2*atan((self.F/self.rhorev(x,y))**(1./self.N))-pi/2
        lon=self.lon_ref+self.theta(x,y)/self.N
        lat=self.convout(lat)
        lon=self.convout(lon)
        return lat,lon
        
    def theta(self,x,y):
        return atan(float(x)/(self.RHO0-y))
        
    def rhorev(self,x,y):
        return self.sgnN*sqrt(x**2+(self.RHO0-y)**2)

    def f(self):
        return cos(self.lat_std1) * tan(pi/4+self.lat_std1/2)**self.N /self.N

    def n(self):
        return log(
               cos(self.lat_std1)/cos(self.lat_std2)
             )/log(
               tan(pi/4+self.lat_std2/2)/tan(pi/4+self.lat_std1/2)
             )

    def rho(self,lat):
        return self.F/tan(pi/4+lat/2)**self.N

    def rho0(self):
        return self.F/tan(pi/4+self.lat_ref/2)**self.N
