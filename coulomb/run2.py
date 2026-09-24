import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import sys, numpy as np, itertools
from multiprocessing import Pool
from kast import logcount
N=int(sys.argv[1]); kind=sys.argv[2]
if kind=="mut":
    white=[(r,c) for r in range(N) for c in range(N) if (r+c)%2==1]
    pairs=list(itertools.combinations(white,2)); rem=[(0,0),(N-1,N-1)]
else:
    white=[(r,c) for r in range(N) for c in range(N) if (r+c)%2==1]
    black=[(r,c) for r in range(N) for c in range(N) if (r+c)%2==0]
    pairs=[(a,b) for a in white for b in black]; rem=[]
def f(ab): return logcount(N,rem,list(ab))
with Pool() as p: L=p.map(f,pairs,chunksize=64)
np.save(f"{kind}{N}.npy",np.array([[a[0],a[1],b[0],b[1],l] for (a,b),l in zip(pairs,L)]))
print("done",kind,N)
