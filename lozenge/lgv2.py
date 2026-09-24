import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
exec(open('pair.py').read().split('x,y=sp.symbols')[0])
import flint
from math import comb
def npath(P,Q):
    dx=Q[0]-P[0]; dy=Q[1]-P[1]
    return comb(dx+dy,dx) if dx>=0 and dy>=0 else 0
def lgvW(a,b,c,k,l):
    src=[(-i,i) for i in range(1,c+1)]+[(a-k,1)]; snk=[(a-j,b+j) for j in range(1,c+1)]+[(a-1,l)]
    M=flint.fmpz_mat([[npath(s,t) for t in snk] for s in src]); return abs(int(M.det()))
def lgvZ(a,b,c):
    comb2=lambda n,k: comb(n,k) if 0<=k<=n else 0
    M=flint.fmpz_mat([[comb2(a+b,b-i+j) for j in range(1,c+1)] for i in range(1,c+1)]); return abs(int(M.det()))
if __name__=='__main__':
    for abc in [(4,4,4),(3,5,4),(5,3,2)]:
        a,b,c=abc; T,Z0=table(a,b,c); assert lgvZ(a,b,c)==Z0
        ok=all(lgvW(a,b,c,k,l)==T[k-1][l-1]*Z0 for k in range(1,a+1) for l in range(1,b+1))
        print(abc,'LGV matches Kasteleyn for all dent pairs:',ok)
