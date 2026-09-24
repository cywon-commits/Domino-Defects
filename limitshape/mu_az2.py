import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np, math, time, sys, json
from var2 import Mesh
from kreg import aztec
res={}
for n in [int(a) for a in sys.argv[1:]]:
    t=time.time(); M=Mesh(aztec(n))
    pure0,h0,exc,r=M.solve(schedule=((50,1e-3),(500,1e-4),(3000,1e-5)),iters=60000)
    v,G,pure,exc=M.obj(h0,1e-5,3000)
    order=M.order; mu=np.array([G[M.idx[p]] for p in order])
    side=[(i,p) for i,p in enumerate(order) if p[0]+p[1] in (n,n+1) and p[0]<=n and p[1]<=n]
    xs=np.array([p[0] for i,p in side]); ms=np.array([mu[i] for i,p in side])
    step={X: ms[xs==X].sum() for X in range(0,n+1) if (xs==X).any()}
    # cumulative from centre X=n/2 outward to X=k : 2*sum mu  vs exact ln[C(n,k)/C(n,n/2)]
    rows=[]
    for k in range(1,n//2):
        lr=2*sum(step.get(X,0) for X in range(k, n//2))
        ex=math.log(math.comb(n,k)/math.comb(n,n//2))
        rows.append((k/n,lr,ex))
    print(f"n={n} S={pure0:.3f} exact {math.log(2)/2*n*n:.3f} ({time.time()-t:.0f}s)")
    for kap,lr,ex in rows[::max(1,len(rows)//8)]: print(f"   kappa={kap:.3f}  2*sum mu={lr:+8.3f}  exact={ex:+8.3f}  ratio={lr/ex if ex else float('nan'):.3f}")
    res[n]=rows
json.dump(res,open('mu_az2_%s.json'%'_'.join(sys.argv[1:]),'w'))
