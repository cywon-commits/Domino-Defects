import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np, math, time, json, sys
from var2 import Mesh
from var_arc import shifted
def house(N,k1,k2): return {(y,x) for y in range(N) for x in range(N) if (x+y)<=2*N-2-k1 and (y-x)<=N-1-k2}
FR=float(sys.argv[1]); N=int(sys.argv[2]); fr=[float(a) for a in sys.argv[3:]]
k=int(round(FR*N)); k-=k%2; cells=house(N,k,k)
t=time.time(); M=Mesh(cells); pure0,h0,exc,r=M.solve(schedule=((50,1e-3),(500,1e-4),(3000,1e-5)),iters=60000)
print("base",N,pure0,"%.0fs"%(time.time()-t),flush=True)
pos={p:i for i,p in enumerate(M.order)}; i0=pos[(0,0)]
for f in fr:
    x=int(round(f*N)); x-=x%2; y=x+N-1-k
    p=(x,y+1) if (x,y+1) in pos else (x,y); i1=pos[p]
    hb=shifted(M,i1,i0,2.0); pure,h,exc,r=M.solve(hb=hb,h0=h0,schedule=((3000,1e-5),),iters=40000)
    print(f"N={N} x/N={x/N:.3f} macro gain {pure-pure0:.3f} ({time.time()-t:.0f}s)",flush=True)
