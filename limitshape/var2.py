import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np, math
from scipy.optimize import minimize
from varsolve import ent_and_grad, boundary_heights
class Mesh:
    """4-triangle-per-cell mesh over a lattice region (cells), boundary heights from lattice (smoothed)."""
    def __init__(s,cells):
        s.cells=set(cells); hb=boundary_heights(s.cells)
        # smooth lattice boundary heights along the boundary (removes flat-edge 0/1 oscillation)
        pts=[];tri=[]; idx={}
        def node(p):
            if p not in idx: idx[p]=len(pts); pts.append(p)
            return idx[p]
        for (y,x) in s.cells:
            a=node((x,y)); b=node((x+1,y)); c=node((x+1,y+1)); d=node((x,y+1)); m=node((x+.5,y+.5))
            tri+= [(a,b,m),(b,c,m),(c,d,m),(d,a,m)]
        s.pts=np.array(pts,float); s.tri=np.array(tri); s.idx=idx
        # boundary nodes = lattice boundary points (centres are always interior for cell-union regions)
        s.bnd=np.array([tuple(p) in hb for p in map(tuple,s.pts)])
        # smoothing: average height of boundary point and its two boundary neighbours (in walk order)
        s.hb_raw=np.array([hb.get(tuple(p),0)/2 for p in map(tuple,s.pts)])
        s.hb=s.hb_raw.copy()
        order=s._walk(hb)
        vals=np.array([hb[p]/2 for p in order]); L=len(order)
        sm=(np.roll(vals,1)+2*vals+np.roll(vals,-1))/4
        for p,v in zip(order,sm): s.hb[idx[p]]=v
        s.order=order
        P=s.pts; A=P[s.tri[:,0]];B=P[s.tri[:,1]];C=P[s.tri[:,2]]
        s.area=0.5*np.abs((B[:,0]-A[:,0])*(C[:,1]-A[:,1])-(C[:,0]-A[:,0])*(B[:,1]-A[:,1]))
        s.Minv=np.linalg.inv(np.stack([B-A,C-A],1))
    def _walk(s,hb):
        # reconstruct boundary order by nearest-neighbour chaining on boundary points (unit steps)
        pts=set(hb); start=min(pts); order=[start]; seen={start}; cur=start
        while True:
            x,y=cur; nxt=None
            for d in ((1,0),(0,1),(-1,0),(0,-1)):
                q=(x+d[0],y+d[1])
                if q in pts and q not in seen and abs(hb[q]-hb[cur])==1: nxt=q; break
            if nxt is None: break
            order.append(nxt); seen.add(nxt); cur=nxt
        return order
    def grads(s,h):
        dA=h[s.tri[:,1]]-h[s.tri[:,0]]; dB=h[s.tri[:,2]]-h[s.tri[:,0]]
        return np.einsum('tij,tj->ti',s.Minv,np.stack([dA,dB],1))
    def obj(s,h,eps,K):
        g=s.grads(h); l1=np.abs(g).sum(1); cap=1-eps
        sc=np.where(l1>cap,cap/np.maximum(l1,1e-300),1.0); gs=g*sc[:,None]
        e,es,et=ent_and_grad(gs[:,0],gs[:,1]); exc=np.maximum(l1-cap,0)
        val=np.sum(s.area*(e/2-K*exc**2))
        cg=s.area[:,None]*(np.stack([es,et],1)/2)-(s.area*2*K*exc)[:,None]*np.sign(g)
        w=np.einsum('tij,ti->tj',s.Minv,cg); G=np.zeros(len(s.pts))
        np.add.at(G,s.tri[:,1],w[:,0]); np.add.at(G,s.tri[:,2],w[:,1]); np.add.at(G,s.tri[:,0],-w[:,0]-w[:,1])
        return val,G,float(np.sum(s.area*e/2)),float(np.sum(s.area*exc))
    def start(s,hb):
        Bp=s.pts[s.bnd]; Bh=hb[s.bnd]; h=np.empty(len(s.pts))
        for i,p in enumerate(s.pts):
            d=np.max(np.abs(Bp-p),axis=1); h[i]=0.5*(np.min(Bh+d)+np.max(Bh-d))
        return h
    def solve(s,hb=None,h0=None,schedule=((50,1e-3),(300,1e-4),(2000,1e-5)),iters=30000):
        hb=s.hb if hb is None else hb
        h=s.start(hb) if h0 is None else h0.copy()
        h[s.bnd]=hb[s.bnd]; free=~s.bnd
        for K,eps in schedule:
            def fg(x):
                hh=h.copy(); hh[free]=x; v,G,_,_=s.obj(hh,eps,K); return -v,-G[free]
            r=minimize(fg,h[free],jac=True,method='L-BFGS-B',options={'maxiter':iters,'maxfun':2*iters,'maxcor':50,'ftol':1e-16,'gtol':1e-10})
            h[free]=r.x
        v,G,pure,exc=s.obj(h,schedule[-1][1],schedule[-1][0])
        return pure,h,exc,r
