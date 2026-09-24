import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np, math
from scipy.special import roots_legendre
class Plane:
    def __init__(s,nq=20000):
        t,w=roots_legendre(nq); th=np.concatenate([(t+1)*np.pi/4,(t+1)*np.pi/4+np.pi/2]); wt=np.concatenate([w,w])*np.pi/4
        c=2*np.cos(th); sg=np.sign(c); sq=np.sqrt(c*c+4)
        s.th,s.wt,s.logz,s.pref=th,wt,np.log(1j*(c-sg*sq)/2),sg/sq
    def T(s,dx,nmax):
        P=(1/np.pi)*s.wt*np.cos(s.th*dx); out=np.empty(nmax+1,complex)
        for a in range(0,nmax+1,48):
            n=np.arange(a,min(a+48,nmax+1))[:,None]; out[a:a+len(n)]=(s.pref[None,:]*np.exp(n*s.logz[None,:]))@P
        return out
def pair_ratio(Pl,d,x=0):
    """black monomer u=(0,x) [x even], white v=(d,x) [d odd]; string = horizontal edges (k,x)-(k,x+1), 0<=k<d"""
    u=(0,x); v=(d,x)
    edges=[]
    for k in range(d):
        a,b=(k,x),(k,x+1); w,bl=(a,b) if (a[0]+a[1])%2 else (b,a); edges.append((w,bl))
    B=[e[1] for e in edges]+[u]; W=[e[0] for e in edges]+[v]
    T={dx:Pl.T(dx,d+2) for dx in (0,1)}
    G=np.array([[T[abs(b[1]-w[1])][abs(b[0]-w[0])] for w in W] for b in B])
    m=d; M=np.zeros((m+1,m+1),complex)
    M[:m,:m]=np.eye(m)-2*G[:m,:m]; M[:m,m]=-2*G[:m,m]; M[m,:m]=G[m,:m]; M[m,m]=G[m,m]
    return abs(np.linalg.det(M)) if m<400 else math.exp(np.linalg.slogdet(M)[1])
Pl=Plane(20000)
from kast import logcount
N=200; L0=logcount(N,[],[])
for d in (1,3,5):
    print("check d=",d, pair_ratio(Pl,d), math.exp(logcount(N,[],[(100,100),(100+d,100)])-L0))
res=[]
for d in [15,31,63,127,255,511,1023,2047]:
    R=pair_ratio(Pl,d); res.append((d,R)); print(d,repr(R), R*d**0.5, flush=True)
import json; json.dump(res,open("bulk_amp.json","w"))
