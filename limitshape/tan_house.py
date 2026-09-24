import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import pickle, numpy as np, math, sys
from tangent_lib import *
from varsolve import ent_and_grad
fn=sys.argv[1]; thr=float(sys.argv[2]) if len(sys.argv)>2 else 0.95
d=pickle.load(open(fn,'rb')); N=d['N']; g=d['g']; C=d['C']
FR=float(fn.split('_')[-1][:-4].rstrip('r')); k=int(round(FR*N)); k-=k%2
l1=np.abs(g).sum(1)
liq=l1<thr; L=C[liq]
PH={(1,0):'E',(-1,0):'W',(0,1):'S',(0,-1):'N'}
def phase_near(P,inward,rad=2.5):
    Q=P+inward*1.5
    sel=np.hypot(C[:,0]-Q[0],C[:,1]-Q[1])<rad
    if not sel.any(): return None,None
    s,t=g[sel].mean(0); Lm=l1[sel].mean()
    if Lm<thr: return 'liquid',g[sel].mean(0)
    return PH.get((int(round(s)),int(round(t)))),None
def mu_at(P,n_out):
    ph,gl=phase_near(P,-n_out)
    if ph is None: return np.nan,ph
    if ph=='liquid':
        e,es,et=ent_and_grad(np.array([gl[0]]),np.array([gl[1]])); F=np.array([es[0],et[0]])
        return float(F@n_out),ph
    t1,t2=tangents(P,L); vals=[]
    for t in (t1,t2):
        F,Q=F_of_dir(t,ph)
        if F is not None: vals.append(float(F@n_out))
    if len(vals)!=1: return (vals[0] if vals else np.nan),ph+('?' if len(vals)!=1 else '')
    return vals[0],ph
c0=N-k   # stair line Y = X + c0 (point coords)
ns=np.array([-1,1])/math.sqrt(2); nl=np.array([-1,0])
# left edge from (0,0) up to junction (0,c0)
ys=np.arange(0.5,c0,1.0); left=[mu_at(np.array([0.0,y]),nl) for y in ys]
xs=np.arange(0.5,k,0.5); stair=[mu_at(np.array([x,x+c0]),ns) for x in xs]
Lsum=np.nansum([m for m,p in left])*1.0
print(fn,'thr',thr,'N',N,'k',k)
print(' left edge: integral mu =',round(Lsum,3),' phases:',''.join((p or '-')[0] for m,p in left))
cum=Lsum; best=(0,-1e9); rows=[]
for x,(m,p) in zip(xs,stair):
    cum+= (m if np.isfinite(m) else 0)*0.5*math.sqrt(2)
    rows.append((x/N,m,p,cum/N))
    if cum/N>best[1]: best=(x/N,cum/N)
for r in rows[::max(1,len(rows)//12)]: print('  x/N=%.3f mu=%+.3f phase=%s  R(x)=%.4f'%(r[0],r[1],r[2],r[3]))
print(' BEST rate R=%.4f at x/N=%.3f'%(best[1],best[0]))
import json
json.dump(dict(rows=[[r[0],r[3]] for r in rows],best=best,left=Lsum/N),open(fn.replace('.pkl','_R.json'),'w'))
