import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np, math, sys, json
from theory import fmap, CORNERS
q={'a1':.75,'h8':.75,'h1':.25,'a8':.25}
def theoryW(N,sites):
    F=[fmap((x+.5)/N,(y+.5)/N) for y,x in sites]
    f=np.array([a for a,b in F]); df=np.array([abs(b) for a,b in F])
    one=df**.25*f.imag**.25
    for c,qq in q.items(): one*=np.abs(f-CORNERS[c])**(-qq)
    P=(np.abs(f[:,None]-f[None,:])*np.abs(f[:,None]-np.conj(f)[None,:]))**0.5
    T=one[:,None]*one[None,:]*P; np.fill_diagonal(T,0); return T
out={}
for N in [int(a) for a in sys.argv[1:]]:
    W=np.load(f"mutW{N}.npy"); S=np.load(f"mutsites{N}.npy"); n=len(S)
    Z=W.sum()/2; Pm=W.sum(1)/(2*Z)        # marginal prob per site (sums to 2 over sites... per defect) 
    T=theoryW(N,[tuple(s) for s in S]); Zt=T.sum()/2; Tm=T.sum(1)/(2*Zt)
    idx={tuple(s):i for i,s in enumerate(S)}; u=1/n
    a1n=Pm[idx[(0,1)]]/u; h1=Pm[idx[(0,N-1)]]/u
    # coarse-grained marginal comparison (blocks of N/8)
    Bk=N//8; Md=np.zeros((8,8)); Mt=np.zeros((8,8))
    for i,(y,x) in enumerate(S): Md[y//Bk,x//Bk]+=Pm[i]; Mt[y//Bk,x//Bk]+=Tm[i]
    r=np.log(Mt/Md)
    # distance distribution
    Y=S[:,0][:,None]-S[:,0][None,:]; X=S[:,1][:,None]-S[:,1][None,:]; D=np.hypot(X,Y)/N
    bins=np.linspace(0,1.42,15)
    hd,_=np.histogram(D,bins,weights=W/(2*Z)); ht,_=np.histogram(D,bins,weights=T/(2*Zt))
    meand=(W*D).sum()/(2*Z); meant=(T*D).sum()/(2*Zt)
    # screening: one defect on each side of anti-diagonal
    side=(S[:,0]+S[:,1])-(N-1); sep=(np.sign(side)[:,None]*np.sign(side)[None,:])<0
    scr=(W*sep).sum()/(2*Z); scrt=(T*sep).sum()/(2*Zt)
    # conditional on v1 at centre
    c=idx[(N//2,N//2+1)]; cd=W[c]/W[c].sum(); ct=T[c]/T[c].sum()
    Cd=np.zeros((8,8)); Ct=np.zeros((8,8))
    for i,(y,x) in enumerate(S): Cd[y//Bk,x//Bk]+=cd[i]; Ct[y//Bk,x//Bk]+=ct[i]
    rc=np.log(Ct/Cd)
    # short-distance repulsion: conditional prob vs distance from centre defect (annuli), compare slope
    dist=np.hypot(S[:,0]-N//2,S[:,1]-(N//2+1))
    ann=[(2,4),(4,8),(8,16)] if N>=48 else [(2,4),(4,8)]
    dens=[]
    for a,b in ann:
        m=(dist>=a)&(dist<b); dens.append(cd[m].sum()/m.sum())
    out[N]=dict(a1=a1n,h1=h1,rms=float(np.sqrt((r**2).mean())),rms_in=float(np.sqrt((r[1:-1,1:-1]**2).mean())),
                meand=meand,meant=meant,scr=scr,scrt=scrt,rms_cond=float(np.sqrt((rc**2).mean())),dens=dens,
                hd=hd.tolist(),ht=ht.tolist(),Md=(Md/Md.mean()).tolist(),Mt=(Mt/Mt.mean()).tolist(),
                Cd=(Cd/Cd.mean()).tolist(),Ct=(Ct/Ct.mean()).tolist())
    print(f"N={N}: a1-adjacent x{a1n:.3f}, h1 corner x{h1:.3f} | marginal rms {out[N]['rms']:.3f} (inner {out[N]['rms_in']:.3f}) | "
          f"<d>/N {meand:.4f} vs {meant:.4f} | screening {scr:.4f} vs {scrt:.4f} | conditional rms {out[N]['rms_cond']:.3f} | "
          f"annulus densities {[round(x*n,3) for x in dens]}")
json.dump({str(k):v for k,v in out.items()},open("mut_an_"+"_".join(sys.argv[1:])+".json","w"))
