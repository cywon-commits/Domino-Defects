import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import pickle, json, numpy as np, math, subprocess
from tangent_lib import *
d=pickle.load(open('fem_az48.pkl','rb')); n=d['N']; g=d['g']; C=d['C']
l1=np.abs(g).sum(1); xh=(C[:,0]-n)/n; yh=(C[:,1]-n)/n
L=np.stack([xh,yh],1)[l1<0.95]; L=L[np.hypot(L[:,0],L[:,1])<0.8]
ks=np.linspace(0.02,0.48,47); mus=[]
for kap in ks:
    P=np.array([-1+kap,-kap]); t1,t2=tangents(P,L); v=None
    for t in (t1,t2):
        F,Q=F_of_dir(t,'W')
        if F is not None: v=-(F[0]+F[1])
    mus.append(abs(v) if v is not None else np.nan)
mus=np.array(mus)
# cumulative from kappa to 1/2
cum=[float(np.trapezoid(np.r_[mus[i:],0.0], np.r_[ks[i:],0.5])) for i in range(len(ks))]
exact_mu=[[float(k),math.log((1-k)/k)] for k in np.linspace(0.01,0.49,97)]
Gc=lambda k: (0.5*math.log(2)) - (-(k*math.log(k)+(1-k)*math.log(1-k)))*1 + 0  # placeholder
cont=[[float(k), math.log(2)-(-(k*math.log(k)+(1-k)*math.log(1-k)))] for k in np.linspace(0.01,0.49,97)]
lr=json.load(open('mu_az2_24_32_48.json'))['48']
naive=[[r[0],-r[1]/48] for r in lr]
out=dict(mu_exact=exact_mu, mu_fem=[[float(k),float(m)] for k,m in zip(ks,mus) if np.isfinite(m)],
         cum_exact=cont, cum_fem=[[float(k),c] for k,c in zip(ks,cum)], cum_naive=naive)
json.dump(out,open('fig_mu.json','w')); print('ok', len(out['mu_fem']))
