import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np, math, json
from rect import Rect, rect_map
res={}
def analyse(H,W,blacks,qmap,tag,extra=None):
    R=Rect(H,W); S=sorted(R.wi); n=len(S)
    Wt=np.array([R.ratio([v],blacks) for v in S]); P=Wt/Wt.sum()
    f,C=rect_map(H,W)
    T=[]
    for (y,x) in S:
        fz,df=f(x+.5,y+.5); v=abs(df)**.25*max(fz.imag,1e-300)**.25
        for key,qq in qmap.items():
            pt=C[key] if isinstance(key,str) else f(key[1],0.0)[0].real
            v*=abs(fz-pt)**(-qq)
        T.append(v)
    T=np.array(T); T/=T.sum()
    nb=8; By=H/nb; Bx=W/nb; Md=np.zeros((nb,nb)); Mt=np.zeros((nb,nb))
    for i,(y,x) in enumerate(S): Md[int(y//By),int(x//Bx)]+=P[i]; Mt[int(y//By),int(x//Bx)]+=T[i]
    r=np.log(Mt/Md); idx={s:i for i,s in enumerate(S)}; u=1/n
    out=dict(rms=float(np.sqrt((r**2).mean())),rms_in=float(np.sqrt((r[1:-1,1:-1]**2).mean())))
    if extra:
        for name,site in extra.items(): out[name]=float(P[idx[site]]/u)
    res[f"{tag}_{H}x{W}"]=out; print(tag,H,W,out,flush=True)
    return P,S
# A: odd x even rectangle, corner (0,0) removed
for H,W in [(15,32),(31,64),(47,96),(63,128)]:
    q={'bl':.75,'br':.25,'tl':-.25,'tr':.25}
    analyse(H,W,[(0,0)],q,"oddeven",extra={'bl_adj':(0,1),'br':(0,W-1),'tl_adj':(H-1,1) if (H-1+1)%2 else (H-2,0),'tr':(H-1,W-1)})
# B: edge mutilation, even square, remove black cell at middle of bottom edge
for N in [32,48,64,96]:
    x0=N//2   # N multiple of 4 -> (0,x0) black
    q={'bl':-.25,'tr':-.25,'br':.25,'tl':.25,(0,x0+0.5):1.0}
    P,S=analyse(N,N,[(0,x0)],q,"edge",extra={'above':(1,x0),'edge1':(0,x0+1),'edge5':(0,x0+5),'edge9':(0,x0+9)})
json.dump(res,open("rect_tests.json","w"))
