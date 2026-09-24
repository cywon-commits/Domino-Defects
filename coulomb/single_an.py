import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np, math, json
from theory import fmap, rho
q={'a1':.75,'h8':-.25,'h1':.25,'a8':.25}
out={}
for N in (32,64,96,128):
    d=np.load(f"single{N}.npy"); L=d[:,2]; P=np.exp(L-L.max()); P/=P.sum()
    S=[(int(y),int(x)) for y,x,_ in d]; n=len(S); u=1/n
    T=np.array([rho(*fmap((x+.5)/N,(y+.5)/N),q) for y,x in S]); T/=T.sum()
    idx={s:i for i,s in enumerate(S)}
    # coarse maps (8x8 blocks of N/8) and (N/16)
    def cg(Pv,nb):
        Bk=N//nb; M=np.zeros((nb,nb))
        for i,(y,x) in enumerate(S): M[y//Bk,x//Bk]+=Pv[i]
        return M
    Md,Mt=cg(P,16),cg(T,16); r=np.log(Mt/Md)
    Y=np.array([s[0] for s in S]); X=np.array([s[1] for s in S])
    dist=np.hypot(X+.5,Y+.5)/N
    quad=(X<N//2)&(Y<N//2)
    out[N]=dict(a1=P[idx[(0,1)]]/u,h1=P[idx[(0,N-1)]]/u,h8=P[idx[(N-1,N-2)]]/u,
        rms=float(np.sqrt((r**2).mean())),rms_in=float(np.sqrt((r[1:-1,1:-1]**2).mean())),
        md=float((P*dist).sum()),mt=float((T*dist).sum()),qd=float(P[quad].sum()),qt=float(T[quad].sum()))
    o=out[N]; print(f"N={N}: rms {o['rms']:.4f} (inner {o['rms_in']:.4f}) | <dist from a1>/N {o['md']:.4f} vs {o['mt']:.4f} | P(a1 quadrant) {o['qd']:.4f} vs {o['qt']:.4f}")
json.dump({str(k):v for k,v in out.items()},open("single_an.json","w"))
