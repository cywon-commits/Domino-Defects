import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np, math, time, sys, json
from var2 import Mesh
def house(N,k1,k2): return {(y,x) for y in range(N) for x in range(N) if (x+y)<=2*N-2-k1 and (y-x)<=N-1-k2}
out={}
for N in [int(a) for a in sys.argv[1:]]:
    t=time.time(); k=N//2; M=Mesh(house(N,k,k))
    pure0,h0,exc,r=M.solve(schedule=((50,1e-3),(500,1e-4),(3000,1e-5)),iters=60000)
    v,G,pure,exc=M.obj(h0,1e-5,3000)
    # pure entropy gradient only
    g=M.grads(h0); 
    from varsolve import ent_and_grad
    pos={p:i for i,p in enumerate(M.order)}; i0=pos[(0,0)]; L=len(M.order)
    mu=np.array([G[M.idx[p]] for p in M.order])
    # walk from a1 backwards (as in house_macro_run: shifted(M,i1,i0) covers i1..i0 cyclic)
    res=[]
    for p,i1 in pos.items():
        x,y=p
        if y-x==N-1-k or y-x==N-k:   # staircase points
            idxs=[(i1+j)%L for j in range((i0-i1)%L+1)]
            gain=2.0*mu[idxs].sum()
            res.append((x/N, gain/N))
    res.sort()
    best=max(res,key=lambda q:q[1])
    print(f"N={N} base S={pure0:.4f} exc={exc:.2e} nit={r.nit} | linear-response max rate={best[1]:.4f} at x/N={best[0]:.3f} ({time.time()-t:.0f}s)",flush=True)
    out[N]=dict(res=res,best=best,mu=mu.tolist(),order=[list(p) for p in M.order])
    json.dump(out,open(f"mu_lr_{'_'.join(sys.argv[1:])}.json","w"))
