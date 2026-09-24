import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import splu
def logcount(N, removed, defects):
    """log(#tilings) of N x N board minus `removed` (boundary cells) minus `defects` (any cells),
    using Kasteleyn weights (horizontal 1, vertical i) with sign flips along a downward dual path
    from every interior removed cell to the bottom boundary."""
    gone=set(removed)|set(defects)
    W=[(r,c) for r in range(N) for c in range(N) if (r+c)%2==1 and (r,c) not in gone]
    B=[(r,c) for r in range(N) for c in range(N) if (r+c)%2==0 and (r,c) not in gone]
    if len(W)!=len(B): return None
    wi={v:i for i,v in enumerate(W)}; bi={v:i for i,v in enumerate(B)}
    flip={}
    for (r,c) in defects:
        if 0<r<N-1 and 0<c<N-1:
            for k in range(r):   # horizontal edges (k,c)-(k,c+1), k<r
                e=((k,c),(k,c+1)); flip[e]=flip.get(e,1)*-1
    rows=[];cols=[];vals=[]
    for (r,c),i in wi.items():
        for dr,dc,w in ((0,1,1),(0,-1,1),(1,0,1j),(-1,0,1j)):
            nb=(r+dr,c+dc)
            if nb in bi:
                e=tuple(sorted([(r,c),nb])) if dr==0 else None
                s=flip.get(e,1) if e else 1
                rows.append(i);cols.append(bi[nb]);vals.append(w*s)
    K=sp.csc_matrix((vals,(rows,cols)),shape=(len(W),len(B)),dtype=complex)
    lu=splu(K)
    return float(np.sum(np.log(np.abs(lu.U.diagonal()))))
