import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np, sys, time
from mutil import Board
N=int(sys.argv[1]); t=time.time()
B=Board(N); Wl=sorted(B.wi)
L=np.array([np.log(B.ratio([v],[(0,0)])) for v in Wl])
np.save(f"single{N}.npy",np.array([[y,x,l] for (y,x),l in zip(Wl,L)])); print("done",N,time.time()-t)
