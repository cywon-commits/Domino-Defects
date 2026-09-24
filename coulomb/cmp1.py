import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np, math
from theory import fmap, rho
from scipy.optimize import least_squares
res={}
for N in (32,64):
    d=np.load(f"single{N}.npy"); L=d[:,2]; P=np.exp(L-L.max()); P/=P.sum()
    res[N]={(int(r),int(c)):p for (r,c,_),p in zip(d,P)}
# corner exponents from data
for N in (8,16,32,64):
    if N in res: P=res[N]
    else:
        W={int(b):(int(hi)<<64)+int(lo) for b,_,hi,lo in (l.split() for l in open(f"O{N}.txt"))}
        t=sum(W.values()); P={divmod(b,N):w/t for b,w in W.items()}
    u=2/N**2
    print(N, "a1-adj %.3f  h1 %.3f  h8-adj %.3f" % (P[(0,1)]/u, P[(0,N-1)]/u, P[(N-1,N-2)]/u))
