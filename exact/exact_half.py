import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import sys, time, json
from exactG import *
from flint import arb, acb, acb_mat, ctx
rs=[int(a) for a in sys.argv[1:]]
Ymax=2*max(rs)+4
cols=plane_table(1,Ymax)
res={}
for r in rs:
    ctx.prec=int(1.4*r*3.33)+300
    pi=arb.pi()
    def G(b,w):
        dx=abs(b[1]-w[1]); return cols[dx][abs(b[0]-w[0])]-cols[dx][b[0]+w[0]+2]
    x=0; v=(r,0); c=(0,0)
    edges=[]
    for k in range(r):
        a,b=(k,x),(k,x+1); ww,bl=(a,b) if (a[0]+a[1])%2 else (b,a); edges.append((ww,bl))
    m=r; t=time.time()
    M=acb_mat(m+1,m+1)
    for j,(wj,bj) in enumerate(edges):
        for l,(wl,bl_) in enumerate(edges):
            val=G(bj,wl).scale(-2).toacb(pi)
            if j==l: val+=1
            M[j,l]=val
        M[j,m]=G(bj,v).scale(-2).toacb(pi)
    for l,(wl,bl_) in enumerate(edges): M[m,l]=G(c,wl).toacb(pi)
    M[m,m]=G(c,v).toacb(pi)
    ctx.prec=256
    d=M.det(); R=abs(d)
    res[r]=R.str(40, radius=False)
    print(r, R.str(30), "time %.1f"%(time.time()-t), flush=True)
json.dump(res,open(f"exact_half_{rs[0]}_{rs[-1]}.json","w"))
