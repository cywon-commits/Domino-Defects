import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import sys, numpy as np
from multiprocessing import Pool
from kast import logcount
N=int(sys.argv[1])
sites=[(r,c) for r in range(N) for c in range(N) if (r+c)%2==1]
def f(v): return logcount(N,[(0,0)],[v])
with Pool() as p: L=p.map(f,sites,chunksize=16)
np.save(f"single{N}.npy",np.array([[r,c,l] for (r,c),l in zip(sites,L)]))
print("done",N)
