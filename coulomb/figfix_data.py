import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
# finer-resolution density maps for Figs. 3 and 4 (adds keys mut96f, odd63s to data/figdata.json)
# needs raw outputs: mutW96.npy, mutsites96.npy (mut_full.py 96) and odd63.npy (odd.py 63) in the working directory
import numpy as np, json
from theory import fmap
from mut_an import theoryW
D=json.load(open('../data/figdata.json'))
# --- mutilated N=96: 2x2 blocks (each block holds exactly 2 white cells) ---
N=96; W=np.load('mutW96.npy'); S=np.load('mutsites96.npy'); n=len(S)
Pm=W.sum(1)/W.sum(); T=theoryW(N,[tuple(s) for s in S]); Tm=T.sum(1)/T.sum()
B=2; nb=N//B; Md=np.zeros((nb,nb)); Mt=np.zeros((nb,nb)); cnt=np.zeros((nb,nb))
for i,(y,x) in enumerate(S): Md[y//B,x//B]+=Pm[i]; Mt[y//B,x//B]+=Tm[i]; cnt[y//B,x//B]+=1
u=1/n; Md=Md/cnt/u; Mt=Mt/cnt/u
print('mut96 2x2: rms log', np.sqrt((np.log(Mt/Md)**2).mean()).round(4), 'max', Md.max().round(2))
D['mut96f']={'data':Md.round(4).tolist(),'theory':Mt.round(4).tolist()}
# --- odd 63x63: per-sublattice maps, no averaging ---
N=63; d=np.load('odd63.npy'); L=d[:,2]; P=np.exp(L-L.max()); P/=P.sum(); u=1/len(P)
pts=[(int(r),int(c)) for r,c,_ in d]
T=np.array([(lambda f,df:(f.imag/abs(df))**0.25)(*fmap((c+.5)/N,(r+.5)/N)) for r,c in pts]); T/=T.sum()
ee=np.full((32,32),np.nan); oo=np.full((31,31),np.nan); tee=ee.copy(); too=oo.copy()
for (r,c),p,t in zip(pts,P,T):
    if r%2==0: ee[r//2,c//2]=p/u; tee[r//2,c//2]=t/u
    else: oo[r//2,c//2]=p/u; too[r//2,c//2]=t/u
# vertical cut through the centre: oo on column 31 (odd rows), ee averaged over columns 30,32 (even rows)
S={k:v/u for k,v in zip(pts,P)}; St={k:v/u for k,v in zip(pts,T)}
prof=[]
for r in range(N):
    if r%2: prof.append([r,S[(r,31)],St[(r,31)]])
    else: prof.append([r,(S[(r,30)]+S[(r,32)])/2,(St[(r,30)]+St[(r,32)])/2])
prof=np.array(prof)
# local sublattice mean and staggering amplitude along the cut
sm=[(r,(prof[r,1]+prof[r+1,1])/2,(prof[r,1]-prof[r+1,1])*(1 if r%2==0 else -1)/(prof[r,1]+prof[r+1,1])) for r in range(31)]
print('sublattice means ee/oo', np.nanmean(ee).round(3), np.nanmean(oo).round(3))
print('cut r, data, theory:'); print(np.round(prof[:12],3))
print('staggering amplitude (ee-oo)/(ee+oo) vs r:', np.round([s[2] for s in sm],3))
rm=np.array([s[1] for s in sm]); tm=np.array([(prof[r,2]+prof[r+1,2])/2 for r in range(31)])
print('pair-mean vs theory ratio:', np.round(rm/tm,3))
D['odd63s']={'ee':np.round(ee,4).tolist(),'oo':np.round(oo,4).tolist(),'tee':np.round(tee,4).tolist(),'too':np.round(too,4).tolist(),'cut':np.round(prof,4).tolist()}
# --- Aztec n=24 single defect (needs az_marg_single24.npy from aztec/az_single.py 24): keep log10(rho/rho_unif) at full precision (cells were rounded to 4 decimals) ---
a=np.load('az_marg_single24.npy')
D['az_single']['log10']=[[int(r),int(c),round(float(np.log10(v)),4)] for r,c,v in a]
json.dump(D,open('../data/figdata.json','w'),allow_nan=True)
