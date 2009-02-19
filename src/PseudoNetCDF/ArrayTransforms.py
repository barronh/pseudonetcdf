"""
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
from PseudoNetCDF.sci_var import PseudoNetCDFVariable

def interior_vertex_func(a,dims=(-1,-2),func=sum):
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

CenterTime=lambda a: interior_vertex_func(a,dims=(0,),func=mean)
CenterLay=lambda a: interior_vertex_func(a,dims=(1,),func=mean)
CenterRow=lambda a: interior_vertex_func(a,dims=(2,),func=mean)
CenterCol=lambda a: interior_vertex_func(a,dims=(3,),func=mean)
CenterRowCol=lambda a: interior_vertex_func(a,dims=(2,3),func=mean)
CenterTimeRowCol=lambda a: interior_vertex_func(a,dims=(0,2,3),func=mean)
CenterCMAQWind=CenterTimeRowCol

def CenterCAMxWind(a):
    out=zeros(a.shape,'f')
    out[:,:,1:-1,1:-1]=CenterRowCol(a[:,:,:-1,:-1])
    out[:,:,[0,-1],1:-1]=CenterCol(a[:,:,[0,-1],:-1])
    out[:,:,1:-1,[0,-1]]=CenterRow(a[:,:,:-1,[0,-1]])
    out[:,:,[0,-1],[0,-1]]=a[:,:,[0,-2],[0,-2]]
    return CenterTime(out)

def CenterCAMxU(a):
    out=zeros(a.shape,'f')
    out[:,:,:,1:-1]=CenterCol(a[:,:,:,:-1])
    return CenterTime(out)

def CenterCAMxV(a):
    out=zeros(a.shape,'f')
    out[:,:,1:-1,:]=CenterRow(a[:,:,:-1,:])
    return CenterTime(out)

def BoundToDiff(a,dim=-1):
    top=[slice(None)]*len(a.shape)
    top[dim]=slice(1,None)
    
    bottom=[slice(None)]*len(a.shape)
    bottom[dim]=slice(None,-1)
    return a[top]-a[bottom]

def CAMxHeightToDepth(a):
    b=a.copy()
    b[:,1:,:,:]=BoundToDiff(a,1)
    return b

def ConvertCAMxTime(date,time,nvars):
    class temp:
        pass
        
        
    f = temp()
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