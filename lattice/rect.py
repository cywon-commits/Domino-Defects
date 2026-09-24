import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np, math, mpmath as mp
class Rect:
    def __init__(s,H,W):
        s.H,s.W=H,W
        Wc=[(y,x) for y in range(H) for x in range(W) if (x+y)%2==1]; Bc=[(y,x) for y in range(H) for x in range(W) if (x+y)%2==0]
        s.wi={v:i for i,v in enumerate(Wc)}; s.bi={v:i for i,v in enumerate(Bc)}
        K=np.zeros((len(Wc),len(Bc)),complex)
        for (y,x),i in s.wi.items():
            for dy,dx,wt in ((0,1,1),(0,-1,1),(1,0,1j),(-1,0,1j)):
                j=s.bi.get((y+dy,x+dx))
                if j is not None: K[i,j]=wt
        s.K=K
        if K.shape[0]==K.shape[1]: s.G=np.linalg.inv(K)
    def string(s,v):
        y,x=v; ks=range(y) if y<s.H/2 else range(y+1,s.H); xl=x if x<s.W-1 else x-1
        return [((k,xl),(k,xl+1)) for k in ks]
    def ratio(s,whites,blacks,black_strings=()):
        cnt={}
        for v in list(whites)+list(black_strings):
            for e in s.string(v): cnt[e]=cnt.get(e,0)^1
        edges=[e for e,c in cnt.items() if c]; ew=[];eb=[]
        for a,b in edges:
            w,bl=(a,b) if (a[0]+a[1])%2==1 else (b,a); ew.append(s.wi[w]); eb.append(s.bi[bl])
        m=len(edges); rows=eb+[s.bi[c] for c in blacks]; cols=ew+[s.wi[v] for v in whites]
        M=s.G[np.ix_(rows,cols)].copy(); M[:m,:]*=-2; M[:m,:m]+=np.eye(m)
        return abs(np.linalg.det(M))
def rect_map(H,W):
    """conformal map of rectangle width W height H (cell units) to UHP: returns function (x,y)->(f,f')"""
    ratio=mp.mpf(2*H)/W
    m=mp.findroot(lambda mm: mp.ellipk(1-mm)/mp.ellipk(mm)-ratio, 0.03 if ratio>1 else 0.5)
    K=mp.ellipk(m); Kp=mp.ellipk(1-m); k=mp.sqrt(m)
    def f(x,y):
        u=-K+2*K*x/W+1j*Kp*y/H
        sn=mp.ellipfun('sn',u,m=m); cn=mp.ellipfun('cn',u,m=m); dn=mp.ellipfun('dn',u,m=m)
        return complex(sn), complex(cn*dn)
    corners={'bl':-1.0,'br':1.0,'tr':float(1/k),'tl':-float(1/k)}
    return f,corners
