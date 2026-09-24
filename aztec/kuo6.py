import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
from kuo2 import setup, strip_forced
for n in range(3,31):
    A,b,SW,NE,I=setup(n); A1,b1,SW1,_,_=setup(n-1); T=A1-SW1
    for cell,sh in (((2*n-2,n-1),None),((n-1,2*n-2),None)):
        R,t=strip_forced(I-{cell})
        found=[(dy,dx) for dy in range(-3,4) for dx in range(-3,4) if {(y+dy,x+dx) for y,x in T}==R]
        if n in (3,4,30): print(n,cell,len(t),found)
        assert found
print('ok')
