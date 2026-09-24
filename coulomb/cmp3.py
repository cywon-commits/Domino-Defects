import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np, sys
from theory import fmap, CORNERS
def run(N,kind,B=4):
    d=np.load(f"{kind}{N}.npy"); L=d[:,4]; P=np.exp(L-L.max()); P/=P.sum()
    A=d[:,:2].astype(int); Bv=d[:,2:4].astype(int)
    sites={tuple(v) for v in np.vstack([A,Bv])}
    F={v:fmap((v[1]+0.5)/N,(v[0]+0.5)/N) for v in sites}
    if kind=="mut": qw={'a1':.75,'h8':.75,'h1':.25,'a8':.25}; ea,eb=-1,-1; qb=qw
    else: qw={'a1':-.25,'h8':-.25,'h1':.25,'a8':.25}; ea,eb=-1,1
    def one(v,e):
        f,df=F[v]; x=abs(df)**.25*f.imag**.25
        for c,q in qw.items(): x*=abs(f-CORNERS[c])**(e*q)
        return x
    T=np.array([one(tuple(a),ea)*one(tuple(b),eb)*(abs(F[tuple(a)][0]-F[tuple(b)][0])*abs(F[tuple(a)][0]-np.conj(F[tuple(b)][0])))**(ea*eb/2) for a,b in zip(A,Bv)])
    T/=T.sum()
    # marginal (coarse-grained, both defects)
    def marg(W):
        M=np.zeros((N//B,N//B))
        for (a,b),w in zip(zip(map(tuple,A),map(tuple,Bv)),W): M[a[0]//B,a[1]//B]+=w; M[b[0]//B,b[1]//B]+=w
        return M
    Md,Mt=marg(P),marg(T); r=np.log(Mt/Md)
    inner=r[1:-1,1:-1]
    dist=np.hypot(*(A-Bv).T)/N
    bins=np.linspace(0,1.42,9)
    hd,_=np.histogram(dist,bins,weights=P); ht,_=np.histogram(dist,bins,weights=T)
    print(f"{kind} N={N}: marginal rms log err {np.sqrt(np.mean(r**2)):.3f} (interior {np.sqrt(np.mean(inner**2)):.3f}); <dist>/N data {np.sum(P*dist):.3f} theory {np.sum(T*dist):.3f}")
    print("   dist histogram data  :",np.round(hd,3)); print("   dist histogram theory:",np.round(ht,3))
    return dict(Md=Md,Mt=Mt,hd=hd,ht=ht,bins=bins)
import pickle
out={}
for N,kind in [(int(a.split(':')[0]),a.split(':')[1]) for a in sys.argv[1:]]:
    out[(kind,N)]=run(N,kind)
pickle.dump(out,open("cmp_pairs_"+"_".join(sys.argv[1:]).replace(':','')+".pkl","wb"))
