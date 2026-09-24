import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np, math
from scipy.optimize import minimize
from scipy.integrate import cumulative_trapezoid
# ---- Lobachevsky L(x) = -int_0^x ln|2 sin t| dt on [0,pi], via fine table ----
import mpmath as _mp
_xs=np.linspace(0,np.pi,40001)
_L=np.array([float(_mp.clsin(2,2*x))/2 for x in _xs])
_Lpi=0.0
def Lob(x): return np.interp(x,_xs,_L) - _Lpi*x/np.pi   # enforce L(pi)=0 (fix discretisation drift)
def ent_and_grad(s,t):
    """CKP domino entropy per domino at slope (s,t), |s|+|t|<=1; returns ent, d/ds, d/dt"""
    lo=np.abs(s); hi=1-np.abs(t)
    for _ in range(55):
        u=(lo+hi)/2
        pE=(u+s)/2; pW=(u-s)/2; pN=(1-u+t)/2; pS=(1-u-t)/2
        f=np.sin(np.pi*pN)*np.sin(np.pi*pS)-np.sin(np.pi*pE)*np.sin(np.pi*pW)
        lo=np.where(f>0,u,lo); hi=np.where(f>0,hi,u)
    u=(lo+hi)/2; pE=(u+s)/2; pW=(u-s)/2; pN=(1-u+t)/2; pS=(1-u-t)/2
    P=[np.clip(p,0,1) for p in (pE,pW,pN,pS)]
    ent=sum(Lob(np.pi*p) for p in P)/np.pi
    eps=1e-300
    ds=-0.5*np.log(np.maximum(np.sin(np.pi*P[0]),eps)/np.maximum(np.sin(np.pi*P[1]),eps))
    dt=-0.5*np.log(np.maximum(np.sin(np.pi*P[2]),eps)/np.maximum(np.sin(np.pi*P[3]),eps))
    return ent,ds,dt
# ---- boundary Thurston heights ----
def boundary_heights(cells):
    """walk region boundary counterclockwise (interior on left); step +1 if the cell on the left is black (x+y even), else -1.
       returns dict lattice point -> Thurston height"""
    edges={}
    for (y,x) in cells:
        # cell occupies [x,x+1]x[y,y+1]; boundary sides oriented with cell on left
        if (y-1,x) not in cells: edges[(x,y)]=((x+1,y),(y,x))       # bottom side, moving +x
        if (y,x+1) not in cells: edges[(x+1,y)]=((x+1,y+1),(y,x)) if (x+1,y) not in edges else edges[(x+1,y)]
        # right side: from (x+1,y) to (x+1,y+1)
    # build explicit list of directed edges (start,end,cell) handling multiple edges per start point
    E=[]
    for (y,x) in cells:
        if (y-1,x) not in cells: E.append(((x,y),(x+1,y),(y,x)))
        if (y,x+1) not in cells: E.append(((x+1,y),(x+1,y+1),(y,x)))
        if (y+1,x) not in cells: E.append(((x+1,y+1),(x,y+1),(y,x)))
        if (y,x-1) not in cells: E.append(((x,y+1),(x,y),(y,x)))
    out={}
    for a,b,c in E: out.setdefault(a,[]).append((b,c))
    start=E[0][0]; h={start:0}; cur=start; prev_dir=None; n=0
    while True:
        opts=out[cur]
        if len(opts)==1: nxt,c=opts[0]
        else:  # pinch point: choose the edge turning most to the right... pick leftmost-turn consistent: prefer continuing
            nxt,c=opts[0]
        step=1 if (c[0]+c[1])%2==0 else -1
        if nxt in h and nxt==start:
            assert h[cur]+step==0, "height not closed: winding %d"%(h[cur]+step)
            break
        h[nxt]=h[cur]+step; cur=nxt; n+=1
        if n>10**7: raise RuntimeError
    return h
class Solver:
    def __init__(s,cells):
        s.cells=set(cells)
        hb=boundary_heights(s.cells)
        pts=set()
        for (y,x) in s.cells: pts|={(x,y),(x+1,y),(x,y+1),(x+1,y+1)}
        s.pts=sorted(pts); idx={p:i for i,p in enumerate(s.pts)}; s.idx=idx
        s.bnd=np.array([p in hb for p in s.pts]); s.hb=np.array([hb.get(p,0)/2 for p in s.pts])   # CKP units
        tri=[]
        for (y,x) in s.cells:
            a,b,c,d=idx[(x,y)],idx[(x+1,y)],idx[(x+1,y+1)],idx[(x,y+1)]
            tri.append((a,b,c)); tri.append((a,c,d))
        s.tri=np.array(tri)
        P=np.array(s.pts,float); s.P=P
        # gradient operator per triangle (area 1/2 each, lattice units)
        A=P[s.tri[:,0]];B=P[s.tri[:,1]];C=P[s.tri[:,2]]
        s.area=0.5*np.abs((B[:,0]-A[:,0])*(C[:,1]-A[:,1])-(C[:,0]-A[:,0])*(B[:,1]-A[:,1]))
        # linear basis gradients
        M=np.stack([B-A,C-A],1)   # (T,2,2) rows = edge vectors
        Minv=np.linalg.inv(M)     # grad = Minv @ [hB-hA, hC-hA]
        s.Minv=Minv
    def feasible_start(s):
        # Lipschitz (l_inf) extension average
        Bp=s.P[s.bnd]; Bh=s.hb[s.bnd]; h=np.empty(len(s.pts))
        for i,p in enumerate(s.P):
            d=np.max(np.abs(Bp-p),axis=1); h[i]=0.5*(np.min(Bh+d)+np.max(Bh-d))
        h[s.bnd]=s.hb[s.bnd]; return h
    def obj(s,hin,hfull,grad=True):
        h=hfull.copy(); h[~s.bnd]=hin
        dA=h[s.tri[:,1]]-h[s.tri[:,0]]; dB=h[s.tri[:,2]]-h[s.tri[:,0]]
        g=np.einsum('tij,tj->ti',s.Minv,np.stack([dA,dB],1))
        sl,tl=g[:,0],g[:,1]
        viol=np.abs(sl)+np.abs(tl)-1
        if np.any(viol>1e-9): 
            return (-1e10*(1+viol[viol>0].sum()), np.zeros(len(s.pts))) if grad else -1e10
        e,es,et=ent_and_grad(sl,tl)
        val=np.sum(s.area*e/2)
        if not grad: return val
        # d val / d g = area/2 * (es,et); dg/dh via Minv
        cg=np.stack([s.area*es/2, s.area*et/2],1)            # (T,2)
        w=np.einsum('tij,ti->tj',s.Minv,cg)                  # coefficients for (dA,dB)
        G=np.zeros(len(s.pts))
        np.add.at(G,s.tri[:,1],w[:,0]); np.add.at(G,s.tri[:,2],w[:,1]); np.add.at(G,s.tri[:,0],-w[:,0]-w[:,1])
        return val,G
    def solve(s,hb_override=None,h0=None,iters=3000):
        hfull=s.hb.copy() if hb_override is None else hb_override.copy()
        if h0 is None: h0=s.feasible_start()
        h0=h0.copy(); h0[s.bnd]=hfull[s.bnd]
        f=lambda x: tuple(-np.array(v) if i==0 else -v[~s.bnd] for i,v in enumerate(s.obj(x,hfull)))
        def fg(x):
            v,G=s.obj(x,hfull); return -v,-G[~s.bnd]
        r=minimize(fg,h0[~s.bnd],jac=True,method='L-BFGS-B',options={'maxiter':iters,'maxcor':30,'ftol':1e-14,'gtol':1e-10})
        hfull[~s.bnd]=r.x; val,G=s.obj(r.x,hfull)
        return val,hfull,G,r

def lip_start(S,hb,alpha):
    Bp=S.P[S.bnd]; Bh=hb[S.bnd]; h=np.empty(len(S.pts))
    for i,p in enumerate(S.P):
        d=alpha*np.max(np.abs(Bp-p),axis=1); h[i]=0.5*(np.min(Bh+d)+np.max(Bh-d))
    h[S.bnd]=hb[S.bnd]; return h
def solve_alpha(S,alpha,hb=None,warm=None,iters=4000):
    hb=(S.hb if hb is None else hb)
    hbA=hb*alpha
    h0=lip_start(S,hbA,alpha)
    if warm is not None:
        # try blending warm start (scaled) with strictly feasible h0
        for beta in (1.0,0.9,0.7,0.5,0.3,0.1,0.0):
            hc=beta*warm+(1-beta)*h0; hc[S.bnd]=hbA[S.bnd]
            if S.obj(hc[~S.bnd],hc,grad=False)>-1e9: h0=hc; break
    return S.solve(hb_override=hbA,h0=h0,iters=iters)

def obj_soft(S,hin,hfull,eps=1e-4,K=200.0):
    h=hfull.copy(); h[~S.bnd]=hin
    dA=h[S.tri[:,1]]-h[S.tri[:,0]]; dB=h[S.tri[:,2]]-h[S.tri[:,0]]
    g=np.einsum('tij,tj->ti',S.Minv,np.stack([dA,dB],1))
    l1=np.abs(g).sum(1); cap=1-eps
    sc=np.where(l1>cap,cap/np.maximum(l1,1e-300),1.0)
    gs=g*sc[:,None]
    e,es,et=ent_and_grad(gs[:,0],gs[:,1])
    exc=np.maximum(l1-cap,0)
    val=np.sum(S.area*(e/2-K*exc**2))
    # gradient: inside -> (es,et)/2 ; outside -> ent gradient projected tangentially (approx: use ent grad at gs) plus penalty
    sg=np.sign(g)
    cg=S.area[:,None]*(np.stack([es,et],1)/2)
    cg-= (S.area*2*K*exc)[:,None]*sg
    w=np.einsum('tij,ti->tj',S.Minv,cg)
    G=np.zeros(len(S.pts))
    np.add.at(G,S.tri[:,1],w[:,0]); np.add.at(G,S.tri[:,2],w[:,1]); np.add.at(G,S.tri[:,0],-w[:,0]-w[:,1])
    return val,G,float(np.sum(S.area*e/2))
def solve_soft(S,hb=None,h0=None,iters=20000,eps=1e-4,K=200.0):
    hfull=(S.hb if hb is None else hb).copy()
    if h0 is None: h0=S.feasible_start()
    h0=h0.copy(); h0[S.bnd]=hfull[S.bnd]
    def fg(x):
        v,G,_=obj_soft(S,x,hfull,eps,K); return -v,-G[~S.bnd]
    r=minimize(fg,h0[~S.bnd],jac=True,method='L-BFGS-B',options={'maxiter':iters,'maxfun':iters*2,'maxcor':50,'ftol':1e-15,'gtol':1e-9})
    hfull[~S.bnd]=r.x; v,G,pure=obj_soft(S,r.x,hfull,eps,K)
    return pure,v,hfull,G,r
