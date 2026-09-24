import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
exec(open('pair.py').read().split('x,y=sp.symbols')[0])
import sympy as sp, itertools
from math import comb
def npath(P,Q):
    dx=Q[0]-P[0]; dy=Q[1]-P[1]
    return comb(dx+dy,dx) if dx>=0 and dy>=0 else 0
def lgv(a,b,c,S,T):
    src=[(-i,i) for i in range(1,c+1)]+[S]; snk=[(a-j,b+j) for j in range(1,c+1)]+[T]
    M=sp.Matrix(c+1,c+1,lambda i,j: npath(src[i],snk[j]))
    return abs(M.det())
a,b,c=3,3,3
T,Z0=table(a,b,c)
Wint=[[T[k][l]*Z0 for l in range(b)] for k in range(a)]
for d1,d2,d3,d4 in itertools.product(range(-2,3),repeat=4):
    ok=True
    for k in range(1,a+1):
        for l in range(1,b+1):
            S=(a-1-k+d1,1+d2); Tt=(a-1+d3,l+d4)
            if lgv(a,b,c,S,Tt)!=Wint[k-1][l-1]: ok=False;break
        if not ok: break
    if ok: print('offsets',d1,d2,d3,d4)
print('done')
