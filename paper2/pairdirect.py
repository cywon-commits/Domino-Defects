"""Direct sparse log-determinant for pair correlations: no cancellation.
logZ(board-{u,v}) = 1/2 log|det K'_{-U}|, K' = K with horizontal string between u and v."""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
from pairsparse import ksparse
def logabsdet(M):
    lu=spl.splu(sp.csc_matrix(M),permc_spec='COLAMD')
    return np.sum(np.log(np.abs(lu.U.diagonal())))+np.sum(np.log(np.abs(lu.L.diagonal())))
class PairD:
    def __init__(s,L,t,r0=None):
        s.L=L; s.t=t; s.r0=L//2 if r0 is None else r0
        s.K=ksparse(L,L,t).tolil(); s.lz=0.5*logabsdet(s.K.tocsc())
    def C(s,c1,c2):
        L=s.L; r0=s.r0; K=s.K.copy()
        for c in range(c1+1,c2):
            i,j=r0*L+c,(r0+1)*L+c; K[i,j]=-K[i,j]; K[j,i]=-K[j,i]
        for c in range(c1,c2):
            if c+1<L:
                i,j=r0*L+c,(r0+1)*L+c+1; K[i,j]=-K[i,j]; K[j,i]=-K[j,i]
        keep=np.setdiff1d(np.arange(L*L),[r0*L+c1,r0*L+c2])
        M=K.tocsc()[keep][:,keep]
        return float(np.exp(0.5*logabsdet(M)-s.lz))
if __name__=='__main__':
    import json,time
    b=json.load(open('../data/paper2_pair_bulk.json'))
    for t in ['0.3','0.5','1.0']:
        T0=time.time(); P=PairD(64,float(t)); ref=dict((int(r),c) for r,c in b[t])
        e=max(abs(P.C(32-(r+1)//2,32-(r+1)//2+r)-ref[r])/ref[r] for r in range(1,33))
        print(t,'max rel err vs dense slogdet: %.1e  (%.1fs)'%(e,time.time()-T0),flush=True)
