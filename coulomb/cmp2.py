import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np, math, pickle
from theory import fmap, rho
from scipy.optimize import least_squares
N=64; B=4
d=np.load(f"single{N}.npy"); L=d[:,2]; P=np.exp(L-L.max()); P/=P.sum()
pts=[(int(r),int(c)) for r,c,_ in d]
F={v:fmap((v[1]+0.5)/N,(v[0]+0.5)/N) for v in pts}
nb=N//B
blocks={}
for v,p in zip(pts,P): blocks.setdefault((v[0]//B,v[1]//B),[]).append((v,p))
keys=sorted(blocks)
data=np.array([sum(p for _,p in blocks[kk]) for kk in keys])
def model(params):
    b,qa1,qh8,qw,D=params
    q={'a1':qa1,'h8':qh8,'h1':qw,'a8':qw}
    rv={v:rho(*F[v],q,b=b,D=D) for v in pts}
    s=sum(rv.values())
    return np.array([sum(rv[v] for v,_ in blocks[kk])/s for kk in keys])
pred_C=model((0.25,0.75,-0.25,0.25,0.25))
pred_N=model((0.25,1.0,0,0,0.25))
pred_B=model((0.0,0.5,0,0,0.25))
def err(pr): 
    r=np.log(pr/data); return np.sqrt(np.mean(r**2)), np.max(np.abs(r))
for name,pr in (("C: corner charges (+3/4,-1/4,+1/4,+1/4)",pred_C),("N: Neumann, a1 only",pred_N),("B: fermion-sqrt",pred_B)):
    print(name, "rms log err %.3f  max %.3f"%err(pr))
fit=least_squares(lambda p: np.log(model(p)/data),x0=[0.25,0.75,-0.25,0.25,0.25])
print("free fit b,q_a1,q_h8,q_w,Delta =",np.round(fit.x,3),"rms %.3f"%np.sqrt(np.mean(fit.fun**2)))
pickle.dump(dict(keys=keys,data=data,C=pred_C,N=pred_N,B=pred_B,nb=nb),open("cmp64.pkl","wb"))
fit2=least_squares(lambda p: np.log(model([p[0],p[1],p[2],p[3],0.25])/data),x0=[0.25,0.75,-0.25,0.25])
print("fit (Delta=1/4 fixed): b,q_a1,q_h8,q_w =",np.round(fit2.x,3),"rms %.3f"%np.sqrt(np.mean(fit2.fun**2)))
# exclude blocks touching the boundary (lattice parity effects)
inner=[i for i,kk in enumerate(keys) if 0<kk[0]<nb-1 and 0<kk[1]<nb-1]
r=np.log(pred_C/data)[inner]; print("model C, interior blocks: rms %.3f max %.3f"%(np.sqrt(np.mean(r**2)),np.abs(r).max()))
r=np.log(pred_C/data); worst=np.argsort(-np.abs(r))[:5]; print("worst blocks",[(keys[i],round(r[i],3)) for i in worst])
