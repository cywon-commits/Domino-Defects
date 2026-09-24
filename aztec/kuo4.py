import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
from kuo2 import setup, Z, math
for n in (4,5,6):
    A,b,SW,NE,I=setup(n); c=(1,n); d=(n,1)
    e=[(n-1+t,2*n-1-t) for t in range(n+1)]
    assert all(v in NE for v in e)
    print('n',n)
    print(' Z(I+b_m+e_t):',[[Z(I|{b[m],e[t]}) for t in range(n+1)] for m in range(n+1)])
    print(' Z(I+e_t-c):',[Z((I|{e[t]})-{c}) for t in range(n+1)])
    print(' Z(I+e_t-d):',[Z((I|{e[t]})-{d}) for t in range(n+1)])
