import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import kuo2
from kuo2 import setup, strip_forced, aztec
for n in range(3,31):
    A,b,SW,NE,I=setup(n); c=(1,n); d=(n,1)
    R,t=strip_forced(I-{c}); A1,b1,SW1,_,_=setup(n-1)
    ok1={(y+1,x) for y,x in A1-SW1}==R
    R,t=strip_forced(I-{d}); ok2={(y,x+1) for y,x in A1-SW1}==R
    R,t=strip_forced(I-{c,d}); ok3={(y+2,x+2) for y,x in aztec(n-2)}==R
    # b_n isolated in I-c
    ok4=all(u not in I-{c} for u in [(1,n),(0,n-1),(0,n+1),(-1,n)])
    assert ok1 and ok2 and ok3 and ok4, n
print('forced reductions verified n=3..30')
