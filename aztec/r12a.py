import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
from kuo_lib import Z, aztec, math
from fractions import Fraction as F
import sys
for n in range(2,8):
    A=aztec(n); ZA=2**(n*(n+1)//2)
    S1=[(n-k,k-1) for k in range(1,n+1)]      # SW_k
    S2=[(n+l-1,l-1) for l in range(1,n+1)]    # from west tip going down-right
    assert all(v in A for v in S1+S2)
    T=[[Z(A-{w,b}) for b in S2] for w in S1]
    tot=sum(map(sum,T))
    print(n, F(tot,ZA), 2**(n-1)-F(1,2))
    if n in (4,5):
        for row in T: print('   ',[F(x,ZA*1)*2**n for x in row])
