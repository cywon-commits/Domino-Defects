import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
from flint import fmpz_mat, fmpq_mat, fmpq
from kreg import aztec
import sys, json
def setup(n):
    A=aztec(n); W=sorted(v for v in A if sum(v)%2); B=sorted(v for v in A if sum(v)%2==0)
    wi={v:i for i,v in enumerate(W)}; bi={v:i for i,v in enumerate(B)}
    M=fmpz_mat(len(W),len(B)); Kval={}
    for (y,x),i in wi.items():
        for dy,dx in ((0,1),(0,-1),(1,0),(-1,0)):
            j=bi.get((y+dy,x+dx))
            if j is not None:
                s=1 if dy==0 else (1 if x%2==0 else -1); M[i,j]=s
    return A,W,B,wi,bi,M
def weights(n):
    A,W,B,wi,bi,M=setup(n); G=M.inv()   # G[b? ] : M is W x B, inv is B x W
    ys=[v[0] for v in A]; ymid=(min(ys)+max(ys))/2
    bnd=lambda v: any((v[0]+a,v[1]+b) not in A for a,b in ((1,0),(-1,0),(0,1),(0,-1)))
    tip=min((v for v in A if sum(v)%2==0), key=lambda v:(v[1],v[0]))
    out={}
    for v in W:
        if bnd(v):
            out[v]=abs(G[bi[tip],wi[v]]); continue
        y,x=v; xr=x+1 if (y,x+1) in A else x-1; xl=min(x,xr)
        step=-1 if y<=ymid else 1; k=y+step; edges=[]
        while (k,xl) in A and (k,xl+1) in A:
            edges.append(((k,xl),(k,xl+1))); k+=step
        ew=[];eb=[]
        for a,b in edges:
            w,bl=(a,b) if (a[0]+a[1])%2==1 else (b,a); ew.append(wi[w]); eb.append(bi[bl])
        m=len(edges); rows=eb+[bi[tip]]; cols=ew+[wi[v]]
        Mx=fmpq_mat(m+1,m+1)
        for r_,rr in enumerate(rows):
            for c_,cc in enumerate(cols):
                g=G[rr,cc]
                Mx[r_,c_]= (fmpq(1 if r_==c_ else 0) - 2*g) if r_<m else g   # K=+1 on horizontal edges
        out[v]=abs(Mx.det())
    return A,out,bnd
if __name__=="__main__":
    res={}
    for n in [int(a) for a in sys.argv[1:]]:
        A,w,bnd=weights(n)
        tot=sum(w.values(),fmpq(0)); interior=sum((x for v,x in w.items() if not bnd(v)),fmpq(0))
        res[n]=dict(total=str(tot),interior=str(interior))
        print(n,"total =",tot,"≈",float(tot),"  interior =",interior,"≈",float(interior),flush=True)
    json.dump(res,open("az_exact_"+"_".join(sys.argv[1:])+".json","w"))
