import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np, math, time, sys
from var2 import Mesh
from kreg import aztec
n=int(sys.argv[1]); sched=eval(sys.argv[2])
t=time.time(); M=Mesh(aztec(n))
pure0,h0,exc,r=M.solve(schedule=sched,iters=60000)
K,eps=sched[-1]
v,G,pure,exc=M.obj(h0,eps,K)
# pure-entropy gradient (no penalty, using capped slopes)
from varsolve import ent_and_grad
g=M.grads(h0); l1=np.abs(g).sum(1); cap=1-eps; sc=np.where(l1>cap,cap/np.maximum(l1,1e-300),1.0); gs=g*sc[:,None]
e,es,et=ent_and_grad(gs[:,0],gs[:,1]); cg=M.area[:,None]*(np.stack([es,et],1)/2)
w=np.einsum('tij,ti->tj',M.Minv,cg); Gp=np.zeros(len(M.pts))
np.add.at(Gp,M.tri[:,1],w[:,0]); np.add.at(Gp,M.tri[:,2],w[:,1]); np.add.at(Gp,M.tri[:,0],-w[:,0]-w[:,1])
order=M.order
def curve(GG):
    mu=np.array([GG[M.idx[p]] for p in order])
    side=[(i,p) for i,p in enumerate(order) if p[0]+p[1] in (n,n+1) and p[0]<=n and p[1]<=n]
    xs=np.array([p[0] for i,p in side]); ms=np.array([mu[i] for i,p in side])
    step={X: ms[xs==X].sum() for X in range(0,n+1)}
    out=[]
    for k in (n//8, n//4, 3*n//8):
        lr=2*sum(step[X] for X in range(k,n//2)); ex=math.log(math.comb(n,k)/math.comb(n,n//2)); out.append(lr/ex)
    return out
print(f"n={n} sched={sched} S={pure0:.3f} exc={exc:.1e} nit={r.nit} | ratio(total G) at k/n=1/8,1/4,3/8: {np.round(curve(G),3)} | ratio(pure) {np.round(curve(Gp),3)} ({time.time()-t:.0f}s)",flush=True)
def local(GG,lab):
    mu=np.array([GG[M.idx[p]] for p in order])
    side=[(i,p) for i,p in enumerate(order) if p[0]+p[1] in (n,n+1) and p[0]<=n and p[1]<=n]
    xs=np.array([p[0] for i,p in side]); ms=np.array([mu[i] for i,p in side])
    step=np.array([ms[xs==X].sum() for X in range(0,n+1)])
    ex=np.array([math.log((n-X)/X) if 0<X<n else np.nan for X in range(0,n+1)])
    W=2; row=[]
    for X in range(2,n//2-1,2):
        a=2*step[X-1:X+2].sum(); b=ex[X-1:X+2].sum(); row.append(f"{X/n:.2f}:{abs(a/b):.2f}")
    print(lab,' '.join(row),flush=True)
local(G,'local total'); local(Gp,'local pure ')
