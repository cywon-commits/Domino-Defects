import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import sys, numpy as np
from multiprocessing import Pool
from kreg import logcount_region, aztec
n=int(sys.argv[1]); kind=sys.argv[2]; A=aztec(n)
if kind=="pair":
    W=sorted(v for v in A if sum(v)%2); B=sorted(v for v in A if sum(v)%2==0)
    jobs=[(a,b) for a in W for b in B]
    def f(ab): return logcount_region(A,list(ab))
else:  # remove west-tip black cell, single white defect
    tip=min((v for v in A if sum(v)%2==0), key=lambda v:(v[1],v[0]))
    A2=A-{tip}; W=sorted(v for v in A2 if sum(v)%2)
    jobs=[(a,tip) for a in W]
    def f(ab): return logcount_region(A2,[ab[0]])
with Pool(2) as p: L=p.map(f,jobs,chunksize=32)
np.save(f"az_{kind}{n}.npy",np.array([[a[0],a[1],b[0],b[1],l] for (a,b),l in zip(jobs,L)]))
print("done az",kind,n)
