import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np, math
def Fcirc(Q):
    x,y=Q
    return np.array([-0.5*math.log(abs((1+2*x)/(1-2*x))), 0.5*math.log(abs((1+2*y)/(1-2*y)))])
ARC={'W':lambda q: q[0]<-abs(q[1]), 'E':lambda q: q[0]>abs(q[1]), 'S':lambda q: q[1]<-abs(q[0]), 'N':lambda q: q[1]>abs(q[0])}
def F_of_dir(d,phase):
    d=np.asarray(d,float); d/=np.linalg.norm(d)
    for sgn in (1,-1):
        Q=sgn*np.array([-d[1],d[0]])/math.sqrt(2)
        if ARC[phase](Q): return Fcirc(Q),Q
    return None,None
def tangents(P,L):
    """L: array of liquid points; returns the two extreme directions (unit vectors) from P."""
    v=L-P; ang=np.arctan2(v[:,1],v[:,0])
    # centre angle
    c=np.angle(np.exp(1j*ang).mean()); rel=np.angle(np.exp(1j*(ang-c)))
    a1=c+rel.min(); a2=c+rel.max()
    return np.array([math.cos(a1),math.sin(a1)]), np.array([math.cos(a2),math.sin(a2)])
