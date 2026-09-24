import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import sys, numpy as np, time
from multiprocessing import Pool
from mutil import Board
N=int(sys.argv[1])
B=Board(N); Wl=sorted(B.wi); n=len(Wl); C=[(0,0),(N-1,N-1)]
def row(i):
    out=np.zeros(n)
    for j in range(i+1,n): out[j]=B.ratio([Wl[i],Wl[j]],C)
    return i,out
if __name__=="__main__":
    t=time.time(); M=np.zeros((n,n))
    order=list(range(n)); order=order[::2]+order[1::2]   # balance
    with Pool(2) as p:
        for k,(i,r) in enumerate(p.imap_unordered(row,order,chunksize=4)):
            M[i]=r
            if k%500==0: print(k,n,"%.0fs"%(time.time()-t),flush=True)
    M=M+M.T
    np.save(f"mutW{N}.npy",M); np.save(f"mutsites{N}.npy",np.array(Wl))
    print("done",N,time.time()-t)
