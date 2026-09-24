import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np, math
from var2 import Mesh
def shifted(M,i0,i1,delta,ramp=3):
    """add delta to boundary heights along walk-order indices i0..i1 (cyclic), linear ramps at both ends"""
    hb=M.hb.copy(); L=len(M.order)
    length=(i1-i0)%L
    for j in range(-ramp,length+ramp+1):
        k=(i0+j)%L
        if j<0: f=(j+ramp)/ramp
        elif j>length: f=(length+ramp-j)/ramp
        else: f=1.0
        f=min(max(f,0),1)
        hb[M.idx[M.order[k]]]+=delta*f
    return hb
def arc_gain(M,h0,pure0,i0,i1,delta):
    hb=shifted(M,i0,i1,delta)
    pure,h,exc,r=M.solve(hb=hb,h0=h0,schedule=((2000,1e-5),))
    return pure-pure0
