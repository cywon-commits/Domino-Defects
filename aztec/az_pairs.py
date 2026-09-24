import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np, sys, time, math, json
from multiprocessing import Pool
from azr import Region
from kreg import aztec
n=int(sys.argv[1]); A=aztec(n); R=Region(A)
Wc=sorted(v for v in A if sum(v)%2); Bc=sorted(v for v in A if sum(v)%2==0)
def row(i):
    w=Wc[i]; return i,np.array([R.ratio([w,b]) for b in Bc])
if __name__=="__main__":
    t=time.time(); M=np.zeros((len(Wc),len(Bc)))
    with Pool(2) as p:
        for i,r in p.imap_unordered(row,range(len(Wc)),chunksize=8): M[i]=r
    np.save(f"azpairW{n}.npy",M); print("done",n,time.time()-t, "total/Z=",M.sum())
