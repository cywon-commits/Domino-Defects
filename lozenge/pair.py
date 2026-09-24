import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
from hexlib import *
from fractions import Fraction as F
import sys, sympy as sp
def table(a,b,c):
    U,D=hexagon(a,b,c); Z0=Z(U,D); A=a+c
    Us=[(A-k,0) for k in range(1,a+1)]; Ds=[(A-1,l-1) for l in range(1,b+1)]
    assert all(u in U for u in Us) and all(d in D for d in Ds)
    return [[F(Z(U-{u},D-{d}),Z0) for d in Ds] for u in Us],Z0
x,y=sp.symbols('x y')
for n in (2,3,4,5,6):
    T,Z0=table(n,n,n)
    print('n',n,'row1',[str(v) for v in T[0]])
    print('     diag',[str(T[i][i]) for i in range(n)], ' sum',sum(map(sum,T)))
