import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np, math, time, sys
from var2 import Mesh
from var_arc import shifted
def house(N,k1,k2): return {(y,x) for y in range(N) for x in range(N) if (x+y)<=2*N-2-k1 and (y-x)<=N-1-k2}
FR=float(sys.argv[1]); N=int(sys.argv[2]); m=int(sys.argv[3]); fr=[float(a) for a in sys.argv[4:]]
Nm=N*m; k=int(round(FR*N)); k-=k%2; km=k*m
t=time.time(); M=Mesh(house(Nm,km,km)); pure0,h0,exc,r=M.solve(schedule=((50,1e-3),(500,1e-4),(3000,1e-5)),iters=60000)
print("base",N,m,pure0,"%.0fs"%(time.time()-t),flush=True)
pos={p:i for i,p in enumerate(M.order)}; i0=pos[(0,0)]
for f in fr:
    x=int(round(f*N)); x-=x%2; xm=x*m; y=xm+Nm-1-km
    p=(xm,y+1) if (xm,y+1) in pos else (xm,y); i1=pos[p]
    hb=shifted(M,i1,i0,2.0*m,ramp=3*m); pure,h,exc,r=M.solve(hb=hb,h0=h0,schedule=((3000,1e-5),),iters=40000)
    print(f"FR={FR} N={N} m={m} x/N={x/N:.3f} gain/m^2 {(pure-pure0)/m**2:.3f} ({time.time()-t:.0f}s)",flush=True)
