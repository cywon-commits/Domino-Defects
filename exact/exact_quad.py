import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import sys, time, json
from exactG import V, ZERO
from flint import fmpq, arb, acb, acb_mat, ctx
def columns_needed(Xset, Ymax):
    Xmax=max(Xset); Ytot=Ymax+Xmax+2
    p={1:fmpq(1,4),-1:fmpq(0)}; q={1:fmpq(0),-1:fmpq(0)}
    for n in range(1,Ytot+2,2):
        p[n+2]=(6*n*p[n]-(n-1)*p[n-2])/(n+1); q[n+2]=(6*n*q[n]-(n-1)*q[n-2]-4)/(n+1)
    col0=[ZERO if n%2==0 else (V(0,-p[n],0,-q[n]) if n%4==1 else V(0,p[n],0,q[n])) for n in range(Ytot+2)]
    g=lambda col,y: col[abs(y)]
    col1=[((V(1) if y==0 else ZERO)-(g(col0,y+1)+g(col0,y-1)).muli()).scale(fmpq(1,2)) for y in range(Ytot+1)]
    keep={}
    if 0 in Xset: keep[0]=col0[:Ymax+1]
    if 1 in Xset: keep[1]=col1[:Ymax+1]
    prev,cur=col0,col1
    for x in range(1,Xmax):
        L=len(cur)-1
        new=[(prev[y]+(g(cur,y+1)+g(cur,y-1)).muli()).neg() for y in range(L)]
        prev,cur=cur,new
        if x+1 in Xset: keep[x+1]=cur[:Ymax+1]
    return keep
res={}
for r in [int(a) for a in sys.argv[1:]]:
    t=time.time()
    x=r+1; v=(r,x); c=(0,0)
    edges=[]
    for k in range(r):
        a,b=(k,x),(k,x+1); ww,bl=(a,b) if (a[0]+a[1])%2 else (b,a); edges.append((ww,bl))
    Bs=[e[1] for e in edges]+[c]; Ws=[e[0] for e in edges]+[v]
    Xset=set(); Ymax=0
    for b in Bs:
        for w in Ws:
            Xset|={abs(b[1]-w[1]), b[1]+w[1]+2}; Ymax=max(Ymax,b[0]+w[0]+2)
    cols=columns_needed(Xset,Ymax)
    def GQ(b,w):
        (yb,xb),(yw,xw)=b,w
        X1,X2=abs(xb-xw),xb+xw+2; Y1,Y2=abs(yb-yw),yb+yw+2
        return cols[X1][Y1]-cols[X2][Y1]-cols[X1][Y2]+cols[X2][Y2]
    t1=time.time()
    ctx.prec=int(5.5*r)+400; pi=arb.pi()
    m=r; M=acb_mat(m+1,m+1)
    for j,bj in enumerate(Bs):
        for l,wl in enumerate(Ws):
            val=GQ(bj,wl)
            if j<m: val=val.scale(-2)
            z=val.toacb(pi)
            if j==l and j<m: z+=1
            M[j,l]=z
    ctx.prec=256
    R=abs(M.det())
    res[r]=R.str(40,radius=False)
    print(r, R.str(30), "table %.1fs total %.1fs"%(t1-t,time.time()-t), flush=True)
    json.dump(res,open("exact_quad.json","w"))
