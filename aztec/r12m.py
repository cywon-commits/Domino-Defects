import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
from kuo_lib import Z, aztec
from math import comb
import time
for n in range(13,21):
    t=time.time(); A=aztec(n); ZA=2**(n*(n+1)//2); ok=True
    for k in range(1,n+1):
        for l in range(1,n+1):
            pred=sum(2**(n*(n+1)//2-1-j)*comb(j,k-1)*comb(j,l-1) for j in range(n))
            if Z(A-{(n-k,k-1),(n+l-1,l-1)})!=pred: ok=False
    # other adjacent pair by symmetry check: S1 x S4 and opposite S1 x S3 zero
    S3=[(2*n-1-i,n+i) for i in range(n)]
    opp=all(Z(A-{(n-1,0),v})==0 for v in S3)
    print(n,ok,opp,round(time.time()-t,1),flush=True)
