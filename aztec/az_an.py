import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np
from kreg import aztec
for n,kind in ((12,"pair"),(24,"single")):
    d=np.load(f"az_{kind}{n}.npy"); L=d[:,4]; P=np.exp(L-L.max()); P/=P.sum()
    A=aztec(n); M={}
    for row,p in zip(d,P):
        a=(int(row[0]),int(row[1])); M[a]=M.get(a,0)+p
        if kind=="pair": b=(int(row[2]),int(row[3])); M[b]=M.get(b,0)+p
    nsite=len(A) if kind=="pair" else len(A)//2
    u=(2 if kind=="pair" else 1)/nsite
    # region classification: normalized coords, arctic circle x^2+y^2 = 1/2 (diamond |x|+|y|<=1)
    def xy(v): return ((v[1]-n+0.5)/n,(v[0]-n+0.5)/n)
    inside=lambda v: sum(t*t for t in xy(v))<0.5
    sites=list(M)
    fin=sum(M[v] for v in sites if inside(v))/(2 if kind=="pair" else 1)
    ain=sum(1 for v in sites if inside(v))/len(sites)
    print(f"Aztec n={n} {kind}: prob mass inside arctic circle {fin:.3f} vs area fraction {ain:.3f}")
    if kind=="pair":
        dist=np.hypot(d[:,0]-d[:,2],d[:,1]-d[:,3]); print("  P(adjacent)=%.3f vs uniform %.3f; <dist>/n=%.3f vs uniform %.3f"%(P[dist==1].sum(),(dist==1).mean(),(P*dist).sum()/n,dist.mean()/n))
        # adjacency split by region
        adj=dist==1
        ins=np.array([inside((int(r[0]),int(r[1]))) for r in d])
        print("  P(adjacent & in frozen)=%.3f  P(adjacent & in liquid)=%.3f"%(P[adj&~ins].sum(),P[adj&ins].sum()))
    # print map (enhancement)
    for r in range(2*n-1,-1,-1):
        print(" ".join(f"{M[(r,c)]/u:4.1f}" if (r,c) in M else ("  . " if (r,c) in A else "    ") for c in range(2*n)))
    np.save(f"az_marg_{kind}{n}.npy",np.array([[v[0],v[1],M[v]/u] for v in sites]))
