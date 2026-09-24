import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np, sys, time, pickle
from var2 import Mesh
def house(N,k1,k2): return {(y,x) for y in range(N) for x in range(N) if (x+y)<=2*N-2-k1 and (y-x)<=N-1-k2}
N=int(sys.argv[1]); FR=float(sys.argv[2]); k=int(round(FR*N)); k-=k%2; cells=house(N,k,k); tag=f'house{N}_{FR}r'
t=time.time(); M=Mesh(cells)
pure,h,exc,r=M.solve(schedule=((50,1e-3),),iters=100000); print('s1',pure,r.nit,r.message,flush=True)
prev=None
for K,eps in ((500,1e-4),(3000,1e-5)):
    for rep in range(6):
        pure,h,exc,r=M.solve(h0=h,schedule=((K,eps),),iters=100000)
        print(K,rep,pure,r.nit,str(r.message)[:40],flush=True)
        if prev is not None and abs(pure-prev)<1e-7*abs(pure): break
        prev=pure
g=M.grads(h); P=M.pts; C=(P[M.tri[:,0]]+P[M.tri[:,1]]+P[M.tri[:,2]])/3
pickle.dump(dict(N=N,cells=cells,g=g,C=C,h=h,pts=P,tri=M.tri,order=M.order,hb=M.hb),open(f'fem_{tag}.pkl','wb'))
print(tag,'done',time.time()-t)
