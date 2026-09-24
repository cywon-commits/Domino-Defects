import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import flint, math, sys
from kreg import aztec
def Z(cells):
    W=sorted(v for v in cells if sum(v)%2); B=sorted(v for v in cells if sum(v)%2==0)
    if len(W)!=len(B): return 0
    iw={v:i for i,v in enumerate(W)}; M=flint.fmpz_mat(len(W),len(B))
    for j,(y,x) in enumerate(B):
        for dy,dx in ((0,1),(0,-1),(1,0),(-1,0)):
            u=(y+dy,x+dx)
            if u in iw: M[iw[u],j]= 1 if dy==0 else (1 if x%2==0 else -1)
    return abs(int(M.det()))
