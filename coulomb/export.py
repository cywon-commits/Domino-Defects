import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np, json, pickle, os
from theory import fmap, rho
out={}
# 1. single defect N=64 : 2x2 block coarse-grained, enhancement over uniform
N=64; B=2
d=np.load("single64.npy"); L=d[:,2]; P=np.exp(L-L.max()); P/=P.sum()
pts=[(int(r),int(c)) for r,c,_ in d]
q={'a1':.75,'h8':-.25,'h1':.25,'a8':.25}
T=np.array([rho(*fmap((c+.5)/N,(r+.5)/N),q) for r,c in pts]); T/=T.sum()
nb=N//B; Md=np.zeros((nb,nb)); Mt=np.zeros((nb,nb))
for (r,c),p,t in zip(pts,P,T): Md[r//B,c//B]+=p; Mt[r//B,c//B]+=t
u=1/(nb*nb)
out['single64']={'data':(Md/u).round(4).tolist(),'theory':(Mt/u).round(4).tolist()}
# 2. corner scaling
cs={}
Ns=[8,16,32,64]+([128] if os.path.exists("single128.npy") else [])
for n in Ns:
    if n in (8,16):
        W={int(b):(int(hi)<<64)+int(lo) for b,_,hi,lo in (l.split() for l in open(f"O{n}.txt"))}
        t=sum(W.values()); Pn={divmod(b,n):w/t for b,w in W.items()}
    else:
        dd=np.load(f"single{n}.npy"); LL=dd[:,2]; PP=np.exp(LL-LL.max()); PP/=PP.sum()
        Pn={(int(r),int(c)):p for (r,c,_),p in zip(dd,PP)}
    uu=2/n**2
    cs[n]=[Pn[(0,1)]/uu, Pn[(0,n-1)]/uu, Pn[(n-1,n-2)]/uu]
out['corner']={str(k):v for k,v in cs.items()}
# 3,4 pairs
for kind in ("mut","reg"):
    import subprocess
    pk=pickle.load(open(f"cmp_pairs_24{kind}.pkl","rb"))[(kind,24)]
    Md,Mt=pk['Md'],pk['Mt']; uu=Md.sum()/Md.size
    out[kind+'24']={'data':(Md/uu).round(4).tolist(),'theory':(Mt/uu).round(4).tolist(),
                    'hd':pk['hd'].round(4).tolist(),'ht':pk['ht'].round(4).tolist(),'bins':pk['bins'].round(3).tolist()}
# 5. boundary vs bulk
out['bnd']=[[1,0.31313,0.24743],[3,0.18687,0.15037],[5,0.11547,0.11128],[9,0.06458,0.08175],[17,0.03213,0.05887],[33,0.01410,0.04114],[65,0.00469,0.02644]]
# 6. aztec
for nm,n in (("pair",12),("single",24)):
    a=np.load(f"az_marg_{nm}{n}.npy")
    out['az_'+nm]={'n':n,'cells':[[int(r),int(c),round(float(v),4)] for r,c,v in a]}
    # include blocked cells of region (for single: other colour + removed tip) as null
from kreg import aztec
out['az_single']['region']=[list(v) for v in aztec(24)]
out['az_pair']['region']=[list(v) for v in aztec(12)]
json.dump(out,open("figdata.json","w"))
print(os.path.getsize("figdata.json"), out['corner'])
