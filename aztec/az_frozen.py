import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np
from kreg import aztec
def edgeprobs(cells):
    W=sorted(v for v in cells if sum(v)%2); B=sorted(v for v in cells if sum(v)%2==0)
    wi={v:i for i,v in enumerate(W)}; bi={v:i for i,v in enumerate(B)}
    K=np.zeros((len(W),len(B)),complex)
    for (r,c),i in wi.items():
        for dr,dc,w in ((0,1,1),(0,-1,1),(1,0,1j),(-1,0,1j)):
            if (r+dr,c+dc) in bi: K[i,bi[(r+dr,c+dc)]]=w
    Ki=np.linalg.inv(K)
    # for each cell: prob its domino goes N,S,E,W
    P={}
    for (r,c),i in wi.items():
        for d,(dr,dc,w) in {'E':(0,1,1),'W':(0,-1,1),'N':(1,0,1j),'S':(-1,0,1j)}.items():
            b=(r+dr,c+dc)
            if b in bi: P[((r,c),d)]=abs(w*Ki[bi[b],i])
    return P
n=40; A=aztec(n)
tip=min((v for v in A if sum(v)%2==0), key=lambda v:(v[1],v[0]))
# white boundary cell in the middle of the SW edge (edge cells: r+c = n-1 region? pick from SW staircase)
sw=[v for v in A if sum(v)%2==1 and v[0]<n and v[1]<n and any((v[0]+a,v[1]+b) not in A for a,b in ((-1,0),(0,-1)))]
sw.sort(key=lambda v:v[1]); mid=sw[len(sw)//2]
print("tip",tip,"defect",mid)
P0=edgeprobs(A); P1=edgeprobs(A-{tip,mid})
def frozen_map(P,cells):
    F={}
    for v in cells:
        if sum(v)%2==1:
            ps=[P.get((v,d),0) for d in 'NSEW']; F[v]=max(ps)
    return F
F0=frozen_map(P0,A); F1=frozen_map(P1,A-{tip,mid})
# liquid = max prob < 0.95 ; compare radius of liquid region along 8 directions in normalized coords
import math
def liquid_radius(F,ang):
    best=0
    for (r,c),f in F.items():
        x=(c-n+0.5)/n; y=(r-n+0.5)/n; a=math.atan2(y,x)
        if abs((a-ang+math.pi)%(2*math.pi)-math.pi)<0.08 and f<0.95: best=max(best,math.hypot(x,y))
    return best
for deg in range(0,360,45):
    a=math.radians(deg); print(f"dir {deg:3d}: liquid extent plain {liquid_radius(F0,a):.3f}  with tip+defect {liquid_radius(F1,a):.3f}  (circle 0.707 on axes/diag scale)")
diff=[(v,F1[v]-F0[v]) for v in F1 if v in F0]
big=[d for d in diff if abs(d[1])>0.05]
print("cells whose frozen-probability changes by >0.05:",len(big),"of",len(diff))
ds=sorted(big,key=lambda d:-abs(d[1]))[:8]; print(ds)
