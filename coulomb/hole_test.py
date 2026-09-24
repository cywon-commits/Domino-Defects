import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np, math, json
from rect import Rect, rect_map
from kast import logcount
# validate interior black hole with string
N=24; R=Rect(N,N); L0=logcount(N,[],[])
for v in [(5,6),(18,3),(12,13)]:
    h=(12,12)
    print("check",v,R.ratio([v],[h],black_strings=[h]), math.exp(logcount(N,[],[h,v])-L0))
res={}
for N in [32,48,64,96]:
    h=(N//2,N//2); R=Rect(N,N); S=sorted(R.wi); n=len(S)
    Wt=np.array([R.ratio([v],[h],black_strings=[h]) for v in S]); P=Wt/Wt.sum()
    f,C=rect_map(N,N); fh,_=f(h[1]+.5,h[0]+.5)
    q={'bl':-.25,'tr':-.25,'br':.25,'tl':.25}
    T=[]
    for (y,x) in S:
        fz,df=f(x+.5,y+.5); v=abs(df)**.25*fz.imag**.25*(abs(fz-fh)*abs(fz-np.conj(fh)))**(-0.5)
        for key,qq in q.items(): v*=abs(fz-C[key])**(-qq)
        T.append(v)
    T=np.array(T); T/=T.sum()
    nb=8; Bk=N//nb; Md=np.zeros((nb,nb)); Mt=np.zeros((nb,nb))
    for i,(y,x) in enumerate(S): Md[y//Bk,x//Bk]+=P[i]; Mt[y//Bk,x//Bk]+=T[i]
    r=np.log(Mt/Md); idx={s:i for i,s in enumerate(S)}; u=1/n
    # hole sits at a block corner for these N; exclude the 4 blocks around it for an 'away' rms
    away=np.ones((nb,nb),bool); away[nb//2-1:nb//2+1,nb//2-1:nb//2+1]=False
    out=dict(rms=float(np.sqrt((r**2).mean())),rms_away=float(np.sqrt((r[away]**2).mean())),
             adj=float(P[idx[(h[0],h[1]+1)]]/u), d5=float(P[idx[(h[0],h[1]+5)]]/u), corner_h1=float(P[idx[(0,N-1)]]/u))
    res[N]=out; print(N,out,flush=True)
json.dump({str(k):v for k,v in res.items()},open("hole_test.json","w"))
