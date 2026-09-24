import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np, sys, time, pickle
from var2 import Mesh
from kreg import aztec
def house(N,k1,k2): return {(y,x) for y in range(N) for x in range(N) if (x+y)<=2*N-2-k1 and (y-x)<=N-1-k2}
kind=sys.argv[1]; N=int(sys.argv[2])
if kind=='az': cells=aztec(N); tag=f'az{N}'
else:
    FR=float(sys.argv[3]); k=int(round(FR*N)); k-=k%2; cells=house(N,k,k); tag=f'house{N}_{FR}'
t=time.time(); M=Mesh(cells); pure0,h0,exc,r=M.solve(schedule=((50,1e-3),(500,1e-4),(3000,1e-5)),iters=60000)
g=M.grads(h0); P=M.pts; C=(P[M.tri[:,0]]+P[M.tri[:,1]]+P[M.tri[:,2]])/3
pickle.dump(dict(N=N,cells=cells,g=g,C=C,h=h0,pts=P,tri=M.tri,order=M.order,hb=M.hb),open(f'fem_{tag}.pkl','wb'))
print(tag,'done',time.time()-t)
