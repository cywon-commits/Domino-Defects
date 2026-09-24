import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import pickle, numpy as np, math
from tangent_lib import *
d=pickle.load(open('fem_az48.pkl','rb')); n=d['N']; g=d['g']; C=d['C']
l1=np.abs(g).sum(1); xh=(C[:,0]-n)/n; yh=(C[:,1]-n)/n
for thr in (0.95,0.99):
    L=np.stack([xh,yh],1)[l1<thr]
    L=L[np.hypot(L[:,0],L[:,1])<0.8]   # drop tip artifacts
    print('thr',thr)
    for kap in (0.05,0.1,0.2,0.3,0.4):
        P=np.array([-1+kap,-kap]); t1,t2=tangents(P,L)
        out=[]
        for t in (t1,t2):
            F,Q=F_of_dir(t,'W')
            out.append(None if F is None else round(-(F[0]+F[1]),4))
        print('  kappa',kap,'mu(t1),mu(t2)=',out,' exact',round(math.log((1-kap)/kap),4))
