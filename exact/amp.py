import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np, math, json, sys
from thmB import Quadrant, ratio
nq=int(sys.argv[1]); Q=Quadrant(nq)
out=[]
for r in [64,96,128,192,256,384,512,768,1024,1536,2048]:
    R=ratio(Q.G,(r,r+1)); rho=math.hypot(r+0.5,r+1.5); out.append([r,rho,R]); print(r,repr(R),flush=True)
json.dump(out,open(f"amp_{nq}.json","w"))
