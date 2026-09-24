import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np, math
from numpy.polynomial.legendre import leggauss
from thmB import ratio
class Half:
    def __init__(s,nq=8000):
        t,w=leggauss(nq); th=np.concatenate([(t+1)*np.pi/4,(t+1)*np.pi/4+np.pi/2]); wt=np.concatenate([w,w])*np.pi/4
        c=2*np.cos(th); sg=np.sign(c); sq=np.sqrt(c*c+4)
        s.th,s.wt,s.zin,s.pref=th,wt,1j*(c-sg*sq)/2,sg/sq
    def G(s,pb,pw):
        yb=np.array([p[0] for p in pb]); xb=np.array([p[1] for p in pb]); yw=np.array([p[0] for p in pw]); xw=np.array([p[1] for p in pw])
        nmax=int(yb.max()+yw.max()+2); F=s.pref[None,:]*s.zin[None,:]**np.arange(nmax+1)[:,None]
        out=np.zeros((len(pb),len(pw)),complex)
        for dx in np.unique(np.abs(xb[:,None]-xw[None,:])):
            T=F@((1/np.pi)*s.wt*np.cos(s.th*dx))
            mask=np.abs(xb[:,None]-xw[None,:])==dx
            D=np.abs(yb[:,None]-yw[None,:]); S=yb[:,None]+yw[None,:]+2
            out[mask]=(T[D]-T[S])[mask]
        return out
H=Half()
from kast import logcount
# check vs big box: boundary cell c=(0,150) removed, defect v=(r,150) r odd (same column), N=300
N=300; L0=logcount(N,[],[])
for r in (1,3,7):
    print("check", r, ratio(H.G,(r,150),corner=(0,150)), math.exp(logcount(N,[(0,150)],[(r,150)])-L0))
prev=None
for r in [3,7,15,31,63,127,255,511]:
    R=ratio(H.G,(r,0),corner=(0,0))
    s=""
    if prev: s="local exp %.5f"%(math.log(R/prev[1])/math.log(r/prev[0]))
    print(r,R,s); prev=(r,R)
# structure: is I+A scalar Toeplitz+Hankel up to gauge?  compare A[j,l] with A[j+2,l+2] and A[j+1,l+1]
y=40; x=0
edges=[]
for k in range(y):
    a,b=(k,x),(k,x+1); w,bl=(a,b) if (a[0]+a[1])%2 else (b,a); edges.append((w,bl))
A=-2*H.G([e[1] for e in edges],[e[0] for e in edges])
Tp=lambda j,l: A[j,l]
print("A[5,8],A[6,9],A[7,10] =",np.round([A[5,8],A[6,9],A[7,10]],6))
print("A[20,3],A[21,4] =",np.round([A[20,3],A[21,4]],6))
