"""Trapping of a forced defect by a fixed charge (square lattice, g = 1/2).
Kasteleyn orientation: horizontal edges east, vertical edges up in even columns, down in odd columns.
Every removed vertex gets a downward string (flip horizontal edges (r',c)-(r',c+1), r' < r).
log Z = 1/2 log|det| from a sparse LU of the string-modified matrix."""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
class Board:
    def __init__(s,L):
        s.L=L; I=[];J=[];V=[];hr=[];hc=[]
        for r in range(L):
            for c in range(L):
                i=r*L+c
                if c+1<L: I+=[i,i+1]; J+=[i+1,i]; V+=[1.,-1.]; hr+=[r,r]; hc+=[c,c]
                if r+1<L:
                    j=i+L
                    if c%2==0: I+=[i,j]; J+=[j,i]
                    else: I+=[j,i]; J+=[i,j]
                    V+=[1.,-1.]; hr+=[-1,-1]; hc+=[-1,-1]
        s.I,s.J,s.V=np.array(I),np.array(J),np.array(V); s.hr,s.hc=np.array(hr),np.array(hc)
    def logZ(s,removed):
        L=s.L; V=s.V.copy()
        for (r,c) in removed: V[(s.hc==c)&(s.hr<r)&(s.hr>=0)]*=-1
        K=sp.csc_matrix((V,(s.I,s.J)),shape=(L*L,L*L))
        keep=np.setdiff1d(np.arange(L*L),[r*L+c for r,c in removed])
        M=K[keep][:,keep].tocsc()
        if (np.diff(M.indptr)==0).any() or (np.bincount(M.indices,minlength=M.shape[0])==0).any(): return -np.inf
        try: lu=spl.splu(M,permc_spec='COLAMD')
        except RuntimeError: return -np.inf
        d=np.abs(lu.U.diagonal())
        return -np.inf if (d==0).any() else 0.5*np.sum(np.log(d))
def brute(L,removed):
    import functools
    cells=frozenset((r,c) for r in range(L) for c in range(L))-frozenset(removed)
    @functools.lru_cache(None)
    def f(S):
        if not S: return 1
        r,c=min(S)
        return sum(f(S-{(r,c),q}) for q in [(r,c+1),(r+1,c)] if q in S)
    return f(cells)
if __name__=='__main__':
    bad=n=0
    for L,fixed in [(6,[(2,3),(1,3),(3,3),(2,2),(2,4),(5,0),(5,2)]),(6,[(0,3),(1,3),(1,2),(1,4),(2,3),(4,1),(4,5)]),(5,[]),(7,[(3,2),(2,2),(4,2),(3,1),(3,3),(6,1),(6,5)])]:
        B=Board(L)
        for v in [(r,c) for r in range(L) for c in range(L) if (r,c) not in fixed]:
            rem=fixed+[v]
            if (L*L-len(rem))%2: continue
            z=brute(L,rem); n+=1
            p=np.exp(B.logZ(rem)) if z>0 else 0
            if abs(p-z)>1e-6*max(z,1): bad+=1
    print('checked',n,'configurations, mismatches',bad)
