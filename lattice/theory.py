import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import mpmath as mp, numpy as np
k=(mp.sqrt(2)-1)**2; m=k**2; K=mp.ellipk(m); Kp=mp.ellipk(1-m)
assert abs(Kp-2*K)<1e-10
def fmap(x,y):
    """unit square (x,y) -> upper half plane; corners (0,0)->-1,(1,0)->1,(1,1)->1/k,(0,1)->-1/k"""
    u=-K+2*K*x+1j*2*K*y
    sn=mp.ellipfun('sn',u,m=m); cn=mp.ellipfun('cn',u,m=m); dn=mp.ellipfun('dn',u,m=m)
    return complex(sn), complex(cn*dn)*float(2*K)   # f, df/dz (z in unit-square coords)
CORNERS={'a1':-1.0,'h1':1.0,'h8':float(1/k),'a8':-float(1/k)}
def rho(fz,dfz,q,b=0.25,D=0.25):
    """Coulomb-gas density for a white defect: |f'|^D (Im f)^b prod |f-f_c|^{-q_c}"""
    v=abs(dfz)**D*max(fz.imag,1e-300)**b
    for c,qc in q.items(): v*=abs(fz-CORNERS[c])**(-qc)
    return v
