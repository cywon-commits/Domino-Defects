import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np, math
from kreg import aztec
class Region:
    def __init__(s,cells):
        s.cells=set(cells)
        W=sorted(v for v in s.cells if sum(v)%2==1); B=sorted(v for v in s.cells if sum(v)%2==0)
        s.wi={v:i for i,v in enumerate(W)}; s.bi={v:i for i,v in enumerate(B)}
        K=np.zeros((len(W),len(B)),complex)
        for (y,x),i in s.wi.items():
            for dy,dx,wt in ((0,1,1),(0,-1,1),(1,0,1j),(-1,0,1j)):
                j=s.bi.get((y+dy,x+dx))
                if j is not None: K[i,j]=wt
        s.G=np.linalg.inv(K)
        ys=[v[0] for v in s.cells]; s.ymid=(min(ys)+max(ys))/2
    def is_boundary(s,v):
        y,x=v; return any((y+a,x+b) not in s.cells for a,b in ((1,0),(-1,0),(0,1),(0,-1)))
    def string(s,v):
        if s.is_boundary(v): return []
        y,x=v; xr=x+1 if (y,x+1) in s.cells else x-1
        xl=min(x,xr); out=[]
        step=-1 if y<=s.ymid else 1
        k=y+step
        while (k,xl) in s.cells and (k,xl+1) in s.cells:
            out.append(((k,xl),(k,xl+1))); k+=step
        return out
    def ratio(s,removed):
        """removed: list of cells (balanced colours). Returns Z(region - removed)/Z(region)."""
        whites=[v for v in removed if sum(v)%2==1]; blacks=[v for v in removed if sum(v)%2==0]
        cnt={}
        for v in removed:
            for e in s.string(v): cnt[e]=cnt.get(e,0)^1
        edges=[e for e,c in cnt.items() if c]; ew=[];eb=[]
        for a,b in edges:
            w,bl=(a,b) if (a[0]+a[1])%2==1 else (b,a); ew.append(s.wi[w]); eb.append(s.bi[bl])
        m=len(edges); rows=eb+[s.bi[c] for c in blacks]; cols=ew+[s.wi[v] for v in whites]
        M=s.G[np.ix_(rows,cols)].copy(); M[:m,:]*=-2; M[:m,:m]+=np.eye(m)
        return abs(np.linalg.det(M))
if __name__=="__main__":
    import subprocess, random
    n=5; A=aztec(n); R=Region(A); N=2*n
    out=[(r,c) for r in range(N) for c in range(N) if (r,c) not in A]
    Z0=float(subprocess.check_output(["./cntf",str(N),str(N)]+[str(t) for v in out for t in v]))
    Wc=[v for v in A if sum(v)%2]; Bc=[v for v in A if sum(v)%2==0]; random.seed(3); bad=0
    for _ in range(60):
        a=random.choice(Wc); b=random.choice(Bc)
        dp=float(subprocess.check_output(["./cntf",str(N),str(N)]+[str(t) for v in out+[a,b] for t in v]))
        r=R.ratio([a,b])
        if abs(r-dp/Z0)>1e-9*max(1,dp/Z0): bad+=1; print(a,b,r,dp/Z0)
    print("mismatches",bad)
