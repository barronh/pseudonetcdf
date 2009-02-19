__doc__ = """
.. _ArrayTransforms
:mod:`utils` -- Array Transforms
================================

.. module:: ArrayTransforms
   :platform: Unix, Windows
   :synopsis: Provides utilities for array shape and type manipulation.
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""
__all__ = ['interior_vertex_func', 'CenterTime', 'CenterLay', 'CenterRow', 'CenterCol', 'CenterRowCol', 'CenterTimeRowCol', 'CenterCMAQWind', 'CenterCAMxWind', 'CenterCAMxU', 'CenterCAMxV', 'BoundToDiff', 'CAMxHeightToDepth', 'ConvertCAMxTime']

HeadURL="$HeadURL$"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy$"
__version__ = RevisionNum

import unittest
from numpy import mean, array,sum,arange,zeros,newaxis
from warnings import warn
from PseudoNetCDF.sci_var import PseudoNetCDFFile, \
                                 PseudoNetCDFVariable

def interior_vertex_func(a,dims=(-1,-2),func=sum):
    """Reduces the array dimensions by 1, by applying
    func to neighboring values.
    
    ex:
        >>> x = array([[1,2,3],
                       [4,5,6],
                       [7,8,9]])
        >>> y = interior_vertex_func(x)
        >>> y
        array([[12, 16],
               [24, 28]])
    
    Note: low index neighbor is the first argument to func
    """
    if func not in (sum,mean):
        warn("Function has only been tested with sum and mean")
    
    if len(dims) not in (1,2,3,4):
        warn("Function has only been tested with up to 4-D")
        
    full_shape=array([slice(None)]*len(a.shape))
    out_array=array(a,copy=True)

    for dim in dims:
        zero_to_minus_one=full_shape.copy()
        one_to_end=full_shape.copy()
        zero_to_minus_one[dim]=slice(0,-1)
        one_to_end[dim]=slice(1,None)
        out_array=func(array([out_array[zero_to_minus_one.tolist()],out_array[one_to_end.tolist()]]),0)
    
    return out_array

def CenterTime(a):
    """
    Center a variable across dimension 0.  By IOAPI and 
    general NetCDF convention, dimension 0 is time.
    """
    return interior_vertex_func(a,dims=(0,),func=mean)
def CenterLay(a):
    """
    Center a variable across dimension 1.  By IOAPI 
    convention, dimension 1 is LAY.
    """
    return interior_vertex_func(a,dims=(1,),func=mean)

def CenterRow(a):
    """
    Center a variable across dimension 2.  By IOAPI 
    convention, dimension 2 is row.
    """
    return interior_vertex_func(a,dims=(2,),func=mean)

def CenterCol(a):
    """
    Center a variable across dimension 3.  By IOAPI
    convention, dimension 3 is col.
    """
    return interior_vertex_func(a,dims=(3,),func=mean)

def CenterRowCol(a):
    """
    Center a variable across dimension 2 and 3.  By IOAPI
    convention, dimension 2 is row and 3 is col.
    """
    return interior_vertex_func(a,dims=(2,3),func=mean)

def CenterTimeRowCol(a):
    """
    Center a variable across dimension 0,2,3.  By IOAPI
    convention, dimension 0 is time, 2 is row, and 3 si col.
    """
    return interior_vertex_func(a,dims=(0,2,3),func=mean)
    
def CenterCMAQWind(a):
    """
    CMAQ wind components are located on the corners of a
    cell.  This function centers those values and returns the 
    time averaged and cell centered value.
    """
    return CenterTimeRowCol(a)

def CenterCAMxWind(a):
    """
    This function needs to be reviewed
    """
    warn("Deprecated")
    out=zeros(a.shape,'f')
    out[:,:,1:-1,1:-1]=CenterRowCol(a[:,:,:-1,:-1])
    out[:,:,[0,-1],1:-1]=CenterCol(a[:,:,[0,-1],:-1])
    out[:,:,1:-1,[0,-1]]=CenterRow(a[:,:,:-1,[0,-1]])
    out[:,:,[0,-1],[0,-1]]=a[:,:,[0,-2],[0,-2]]
    return CenterTime(out)

def CenterCAMxU(a):
    """
    Interpolate staggered U component of CAMx wind
    to cell center
    """
    out=zeros(a.shape,'f')
    out[:,:,:,1:-1]=CenterCol(a[:,:,:,:-1])
    return CenterTime(out)

def CenterCAMxV(a):
    """
    Interpolate staggered V component of CAMx wind
    to cell center
    """
    out=zeros(a.shape,'f')
    out[:,:,1:-1,:]=CenterRow(a[:,:,:-1,:])
    return CenterTime(out)

def BoundToDiff(a,dim=-1):
    """
    Calculate the first order difference between
    neighbors along dimension
    """
    from numpy import diff
    
    return diff(zf, axis = dim)

def CAMxHeightToDepth(a):
    """
    Use full height array to calculate the depth
    of each model layer
    """
    b=a.copy()
    b[:,1:,:,:]=BoundToDiff(a,1)
    return b

def ConvertCAMxTime(date,time,nvars=1):
    """
    Use camx date and time arrays to produce an
    IOAPI standard TFLAG variable
    """
    f = PseudoNetCDFFile()
    f.dimensions = {'TSTEP': date.shape[0], 'VAR': nvars, 'DATE-TIME': 2}
    
    a=array([date,time],dtype='i').swapaxes(0,1)
    if len(a.shape)==2:
        a=a[:,newaxis,:]
    date=a[:,:,0]
    if (date<70000).any():
        date+=2000000
    else:
        date+=1900000
    time=a[:,:,1]
    while not (time==0).all() and time.max()<10000:
        time*=100
    a=PseudoNetCDFVariable(f,'TFLAG','i',('TSTEP','VAR','DATE-TIME'),values=a[:,[0],:].repeat(nvars,1))
    a.units='DATE-TIME'.ljust(16)
    a.long_name='TFLAG'.ljust(16)
    a.var_desc=a.long_name
    return a
    
            
    
class TestInteriorVertex(unittest.TestCase):
    def setUp(self):
        self.A2D=arange(1,17).reshape(4,4)
        self.A3D=arange(1,17).reshape(1,4,4).repeat(4,0)
        self.A4D=arange(1,17).reshape(1,1,4,4).repeat(4,0).repeat(4,1)
        self.mean2D=array([[  3.5,   4.5,   5.5],
            [  7.5,   8.5,   9.5],
            [ 11.5,  12.5,  13.5]])
        
        self.sum2D=array([[14, 18, 22],
            [30, 34, 38],
            [46, 50, 54]])

    def runTest(self):
        pass
        
    def test2D2DSum(self):
        self.assert_((interior_vertex_func(self.A2D,dims=(-1,-2),func=sum)==self.sum2D).all())

    def test2D2DMean(self):
        self.assert_((interior_vertex_func(self.A2D,dims=(-1,-2),func=mean)==self.mean2D).all())

    def test3D2DSum(self):
        self.assert_((interior_vertex_func(self.A3D,dims=(-1,-2),func=sum)==self.sum2D.reshape(1,3,3)).all())

    def test3D2DMean(self):
        self.assert_((interior_vertex_func(self.A3D,dims=(-1,-2),func=mean)==self.mean2D.reshape(1,3,3)).all())

    def test3D3DSum(self):
        self.assert_((interior_vertex_func(self.A3D,dims=(-1,-2,-3),func=sum)==(self.sum2D.reshape(1,3,3)*2)).all())

    def test3D3DMean(self):
        self.assert_((interior_vertex_func(self.A3D,dims=(-1,-2,-3),func=mean)==self.mean2D.reshape(1,3,3)).all())

    def test4D2DSum(self):
        self.assert_((interior_vertex_func(self.A4D,dims=(-1,-2),func=sum)==self.sum2D.reshape(1,1,3,3)).all())

    def test4D2DMean(self):
        self.assert_((interior_vertex_func(self.A4D,dims=(-1,-2),func=mean)==self.mean2D.reshape(1,1,3,3)).all())

    def test4D3DSum(self):
        self.assert_((interior_vertex_func(self.A4D,dims=(-1,-2,-3),func=sum)==(self.sum2D.reshape(1,1,3,3)*2)).all())

    def test4D3DMean(self):
        self.assert_((interior_vertex_func(self.A4D,dims=(-1,-2,-3),func=mean)==self.mean2D.reshape(1,1,3,3)).all())

    def test4D4DSum(self):
        self.assert_((interior_vertex_func(self.A4D,dims=(-1,-2,-3,-4),func=sum)==(self.sum2D.reshape(1,1,3,3)*4)).all())

    def test4D4DMean(self):
        self.assert_((interior_vertex_func(self.A4D,dims=(-1,-2,-3,-4),func=mean)==self.mean2D.reshape(1,1,3,3)).all())


if __name__=='__main__':
    unittest.main()