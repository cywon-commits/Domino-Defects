import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import math,itertools
from kast import logcount
N=8
# single defect check vs DP (O8.txt)
dp={int(b):int(hi)*2**64+int(lo) for b,_,hi,lo in (l.split() for l in open("O8.txt"))}
bad=0
for b,w in dp.items():
    r,c=divmod(b,N); lc=logcount(N,[(0,0)],[(r,c)])
    if abs(lc-math.log(w))>1e-6: bad+=1
print("single-defect mismatches:",bad,"of",len(dp))
# mutilated pairs check vs W8
bad=0;n=0
for l in open("W8.txt"):
    a,b,w=l.split(); a,b=int(a),int(b); w=float(w)
    lc=logcount(N,[(0,0),(7,7)],[divmod(a,N),divmod(b,N)]); n+=1
    if abs(lc-math.log(w))>1e-6: bad+=1
print("mutilated mismatches:",bad,"of",n)
bad=0;n=0
for l in open("F8.txt"):
    a,b,w=l.split(); a,b=int(a),int(b); w=float(w)
    lc=logcount(N,[],[divmod(a,N),divmod(b,N)]); n+=1
    if abs(lc-math.log(w))>1e-6: bad+=1
print("regular-board mismatches:",bad,"of",n)
