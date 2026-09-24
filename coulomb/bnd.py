import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import math
from kast import logcount
N=160; L0=logcount(N,[],[])
prev=None
print("d   P_boundary   slope   | P_bulk   slope")
for d in [1,3,5,9,17,33,65]:
    c1=(N-d)//2; c2=c1+d
    if (0+c1)%2==(0+c2)%2: c2+=1; d+=1
    lb=logcount(N,[],[(0,c1),(0,c2)])-L0
    r=N//2; lm=logcount(N,[],[(r,c1),(r,c2)])-L0
    s=""
    if prev: s=f"{(lb-prev[1])/math.log(d/prev[0]):7.3f}      {(lm-prev[2])/math.log(d/prev[0]):7.3f}"
    print(d, f"{math.exp(lb):.5f}", f"{math.exp(lm):.5f}", s); prev=(d,lb,lm)
