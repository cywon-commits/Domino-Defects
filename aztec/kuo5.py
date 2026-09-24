import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
from kuo2 import setup, strip_forced, aztec, Z, math
for n in range(3,31):
    A,b,SW,NE,I=setup(n); c=(1,n); d=(n,1)
    e=[(n-1+t,2*n-1-t) for t in range(n+1)]
    sig=lambda v:(2*n-1-v[1],2*n-1-v[0])
    assert {sig(v) for v in I}==I and {sig(v) for v in b}==set(e)
    assert [sig(v) for v in b]==e[::-1] or [sig(v) for v in b]==e
    if n<=7:
        K=2**((n-1)*(n-2)//2)
        for t in range(n+1):
            assert Z((I|{e[t]})-{c})==(K*math.comb(n-1,t-1) if t>0 else 0)
            assert Z((I|{e[t]})-{d})==K*math.comb(n-1,t)
        for m in range(n+1):
            assert sum(Z(I|{b[m],e[t]}) for t in range(n+1))==2**(n*(n-1)//2)*math.comb(n,m)
print([ (v,) for v in []], 'sig check ok; e-terms ok n<=7; NE switch sum ok')
n=5;A,b,SW,NE,I=setup(n);sig=lambda v:(2*n-1-v[1],2*n-1-v[0]);print([sig(v) for v in b]);print([(n-1+t,2*n-1-t) for t in range(n+1)])
