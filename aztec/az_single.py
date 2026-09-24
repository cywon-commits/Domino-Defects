import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np, json, time, sys
from azr import Region
from kreg import aztec
res={}
for n in [int(a) for a in sys.argv[1:]]:
    t=time.time(); A=aztec(n); R=Region(A)
    tip=min((v for v in A if sum(v)%2==0), key=lambda v:(v[1],v[0]))
    Wc=sorted(v for v in A if sum(v)%2==1)
    w=np.array([R.ratio([tip,v]) for v in Wc]); P=w/w.sum()
    xy=np.array([((c-n+0.5)/n,(r-n+0.5)/n) for r,c in Wc])
    liquid=(xy**2).sum(1)<0.5
    bnd=np.array([R.is_boundary(v) for v in Wc])
    sw=np.array([R.is_boundary(v) and xy[i,0]<0 and xy[i,1]<0 for i,v in enumerate(Wc)])  # SW quadrant boundary
    res[n]=dict(Pb=float(P[bnd].sum()),Psw=float(P[sw].sum()),Pliq=float(P[liquid].sum()),Pfrozen_int=float(P[~bnd&~liquid].sum()),
                Zratio=float(w.sum()),
                sw_profile=[[float(-xy[i,1]-(-0.5/n)), float(P[i])] for i in np.where(sw)[0]])
    np.save(f"azs{n}.npy",np.column_stack([np.array(Wc),P]))
    print(n, {k:round(v,5) for k,v in res[n].items() if k!='sw_profile'}, "%.0fs"%(time.time()-t), flush=True)
json.dump({str(k):v for k,v in res.items()},open("az_single_"+"_".join(sys.argv[1:])+".json","w"))
