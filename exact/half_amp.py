import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np, math, json
from scipy.special import roots_legendre
from thmB import ratio
class Half:
    def __init__(s,nq=20000):
        t,w=roots_legendre(nq); th=np.concatenate([(t+1)*np.pi/4,(t+1)*np.pi/4+np.pi/2]); wt=np.concatenate([w,w])*np.pi/4
        c=2*np.cos(th); sg=np.sign(c); sq=np.sqrt(c*c+4)
        s.th,s.wt,s.logz,s.pref=th,wt,np.log(1j*(c-sg*sq)/2),sg/sq
    def T(s,dx,nmax):
        P=(1/np.pi)*s.wt*np.cos(s.th*dx); out=np.empty(nmax+1,complex)
        for a in range(0,nmax+1,48):
            n=np.arange(a,min(a+48,nmax+1))[:,None]; out[a:a+len(n)]=(s.pref[None,:]*np.exp(n*s.logz[None,:]))@P
        return out
    def G(s,pb,pw):
        yb=np.array([p[0] for p in pb]); xb=np.array([p[1] for p in pb]); yw=np.array([p[0] for p in pw]); xw=np.array([p[1] for p in pw])
        nmax=int(yb.max()+yw.max()+2); out=np.zeros((len(pb),len(pw)),complex)
        DX=np.abs(xb[:,None]-xw[None,:]); D=np.abs(yb[:,None]-yw[None,:]); S=yb[:,None]+yw[None,:]+2
        for dx in np.unique(DX):
            T=s.T(dx,nmax); mask=DX==dx; out[mask]=(T[D]-T[S])[mask]
        return out
H=Half(20000); res=[]
for r in [63,95,127,191,255,383,511,767,1023,1535,2047]:
    R=ratio(H.G,(r,0),corner=(0,0)); res.append([r,R]); print(r,repr(R),flush=True)
json.dump(res,open("half_amp.json","w"))
