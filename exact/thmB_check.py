import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np, math
from thmB import Quadrant, ratio
from kast import logcount
import scipy.sparse as sp
# 1) verify determinant formula on a finite box using dense inverse of the full operator M
N=20
idx={(y,x):i for i,(y,x) in enumerate((y,x) for y in range(N) for x in range(N))}
M=np.zeros((N*N,N*N),complex)
for (y,x),i in idx.items():
    for dy,dx,wgt in ((0,1,1),(0,-1,1),(1,0,1j),(-1,0,1j)):
        j=idx.get((y+dy,x+dx))
        if j is not None: M[i,j]=wgt
Mi=np.linalg.inv(M)
Gbox=lambda Bs,Ws: Mi[np.ix_([idx[b] for b in Bs],[idx[w] for w in Ws])]
L0=logcount(N,[],[])
for v in [(3,4),(5,6),(4,7),(7,8)]:
    print("box check",v, ratio(Gbox,v), math.exp(logcount(N,[(0,0)],[v])-L0))
# 2) quadrant vs large box
Q=Quadrant(3000)
L0=logcount(300,[],[])
for v in [(0,1),(0,3),(3,4),(6,7),(10,11)]:
    print("quadrant",v, ratio(Q.G,v), "box300", math.exp(logcount(300,[(0,0)],[v])-L0))
