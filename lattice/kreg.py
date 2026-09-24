import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import splu
def logcount_region(cells, defects):
    """cells: set of (r,c) forming a simply connected region; defects removed with dual-path sign flips."""
    gone=set(defects)
    W=[v for v in cells if (v[0]+v[1])%2==1 and v not in gone]
    B=[v for v in cells if (v[0]+v[1])%2==0 and v not in gone]
    if len(W)!=len(B): return None
    wi={v:i for i,v in enumerate(W)}; bi={v:i for i,v in enumerate(B)}
    flip={}
    for (r,c) in defects:
        for k in range(r):
            e=((k,c),(k,c+1)); flip[e]=flip.get(e,1)*-1
    rows=[];cols=[];vals=[]
    for (r,c),i in wi.items():
        for dr,dc,w in ((0,1,1),(0,-1,1),(1,0,1j),(-1,0,1j)):
            nb=(r+dr,c+dc)
            if nb in bi:
                s=flip.get(tuple(sorted([(r,c),nb])),1) if dr==0 else 1
                rows.append(i);cols.append(bi[nb]);vals.append(w*s)
    K=sp.csc_matrix((vals,(rows,cols)),shape=(len(W),len(B)),dtype=complex)
    lu=splu(K)
    return float(np.sum(np.log(np.abs(lu.U.diagonal()))))
def aztec(n):
    return {(r,c) for r in range(2*n) for c in range(2*n) if abs(c-n+0.5)+abs(r-n+0.5)<=n}
