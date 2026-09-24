import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
from hexlib import *
from fractions import Fraction as F
from math import comb
import itertools
def row1(a,b,c):
    U,D=hexagon(a,b,c); Z0=Z(U,D); A=a+c
    u=(A-1,0); Ds=[(A-1,l-1) for l in range(1,b+1)]
    vals=[F(1)]+[F(Z(U-{u},D-{d}),Z0) for d in Ds]+[F(0)]
    return [vals[i]-vals[i+1] for i in range(b+1)]
for abc in [(3,4,5),(4,3,5),(5,4,3),(2,5,3),(4,4,2)]:
    w=row1(*abc); a,b,c=abc
    # guess: w_m ∝ C(p+m,m) C(q+b-m, b-m)
    found=[]
    for p in range(0,12):
        for q in range(0,12):
            g=[F(comb(p+m,m)*comb(q+b-m,b-m)) for m in range(b+1)]
            s=sum(g); 
            if all(F(x)==y/s for x,y in zip(w,g)): found.append((p,q))
    print(abc,[str(x) for x in w],'fits C(p+m,m)C(q+b-m,b-m) with (p,q)=',found)
