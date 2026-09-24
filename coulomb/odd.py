import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import sys, numpy as np, math, subprocess
from multiprocessing import Pool
from kast import logcount
def run(N):
    sites=[(r,c) for r in range(N) for c in range(N) if (r+c)%2==0]
    with Pool(2) as p: L=p.starmap(logcount,[(N,[],[v]) for v in sites],chunksize=16)
    return sites,np.array(L)
if __name__=="__main__":
    # verify vs DP
    N=7; s,L=run(N); bad=0
    for v,l in zip(s,L):
        dp=float(subprocess.check_output(["./cntf",str(N),str(N),str(v[0]),str(v[1])]))
        if abs(math.log(dp)-l)>1e-6: bad+=1
    print("7x7 mismatches",bad)
    for N in map(int,sys.argv[1:]):
        s,L=run(N); np.save(f"odd{N}.npy",np.array([[r,c,l] for (r,c),l in zip(s,L)])); print("done",N,flush=True)
