import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np, math
from theory import fmap
for N in (7,15,31,63):
    d=np.load(f"odd{N}.npy"); L=d[:,2]; P=np.exp(L-L.max()); P/=P.sum()
    S={(int(r),int(c)):p for (r,c,_),p in zip(d,P)}; u=1/len(S); m=N//2
    ee=sum(p for (r,c),p in S.items() if r%2==0)/sum(1 for (r,c) in S if r%2==0)
    oo=sum(p for (r,c),p in S.items() if r%2==1)/sum(1 for (r,c) in S if r%2==1)
    print(f"N={N}: corner {S[(0,0)]/u:.3f}  edge-mid {S[(0,m if m%2==0 else m-1)]/u:.3f}  centre {S[(m,m) if (2*m)%2==0 else (m,m)]/u:.3f}  max {max(S.values())/u:.3f}  min {min(S.values())/u:.3f} | mean on (even,even) {ee/u:.3f} vs (odd,odd) {oo/u:.3f}")
    if N==7:
        for r in range(N-1,-1,-1): print("   "+" ".join(f"{S[(r,c)]/u:5.2f}" if (r,c) in S else "   . " for c in range(N)))
# theory comparison N=63, 3x3 blocks... use (even,even)+(odd,odd) together, blocks of 3? use 7x7 coarse (9 blocks of 7)
N=63; d=np.load("odd63.npy"); L=d[:,2]; P=np.exp(L-L.max()); P/=P.sum()
pts=[(int(r),int(c)) for r,c,_ in d]
T=[]
for r,c in pts:
    f,df=fmap((c+.5)/N,(r+.5)/N); T.append((f.imag/abs(df))**0.25)
T=np.array(T); T/=T.sum()
B=7; nb=N//B; Md=np.zeros((nb,nb)); Mt=np.zeros((nb,nb))
for (r,c),p,t in zip(pts,P,T): Md[r//B,c//B]+=p; Mt[r//B,c//B]+=t
rr=np.log(Mt/Md); print("theory rho ∝ (Im f/|f'|)^(1/4): rms log err %.3f, interior %.3f, max %.3f"%(np.sqrt((rr**2).mean()),np.sqrt((rr[1:-1,1:-1]**2).mean()),np.abs(rr).max()))
u=1/nb**2; print("coarse data (x uniform):");print(np.round(Md[::-1]/u,2)); print("theory:");print(np.round(Mt[::-1]/u,2))
# edge exponent: average over parity along the middle column, distance from bottom edge
col=[ (r,S) for r in range(N)]
S={(r,c):p for (r,c),p in zip(pts,P)}
m=N//2
prof=[(S.get((r,m),0)+S.get((r,m-1),0)+S.get((r,m+1),0)) for r in range(N)]
print("middle-column profile rows 0..10 (x centre):",np.round(np.array(prof[:11])/prof[m],3))
