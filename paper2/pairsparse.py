"""Pair correlation via sparse LU (same algebra as pairfast.py, no dense K).
C(r) = Z(board-{u,v})/Z(board), u=(r0,c1), v=(r0,c2), horizontal string between them."""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
def ksparse(Lx, Ly, t):
    I=[];J=[];V=[]
    up=lambda c: 1 if c%2==0 else -1
    def add(a,b,w): I.extend([a,b]); J.extend([b,a]); V.extend([w,-w])
    for r in range(Ly):
        for c in range(Lx):
            i=r*Lx+c
            if c+1<Lx: add(i,i+1,1.0)
            if r+1<Ly:
                j=i+Lx
                if up(c)==1: add(i,j,1.0)
                else: add(j,i,1.0)
            if t!=0 and c+1<Lx and r+1<Ly:
                j=i+Lx+1
                if up(c+1)==1: add(i,j,t)
                else: add(j,i,t)
    n=Lx*Ly; return sp.csc_matrix((V,(I,J)),shape=(n,n))
class Pair:
    def __init__(s,L,t,r0=None):
        s.L=L; s.t=t; s.r0=L//2 if r0 is None else r0
        K=ksparse(L,L,t); S=np.r_[s.r0*L:(s.r0+2)*L]; s.S=S
        s.KS=K[S][:,S].toarray()
        lu=spl.splu(K,permc_spec='COLAMD')
        E=np.zeros((L*L,len(S))); E[S,np.arange(len(S))]=1
        s.G=lu.solve(E)[S,:]; del lu
    def C(s,c1,c2):
        L=s.L; loc=lambda r,c:(r-s.r0)*L+c
        D=np.zeros_like(s.KS)
        for c in range(c1,c2+1):
            if c1<c<c2:
                i,j=loc(s.r0,c),loc(s.r0+1,c); D[i,j]-=2*s.KS[i,j]; D[j,i]-=2*s.KS[j,i]
            if c<c2 and c+1<L:
                i,j=loc(s.r0,c),loc(s.r0+1,c+1); D[i,j]-=2*s.KS[i,j]; D[j,i]-=2*s.KS[j,i]
        idx=np.where(np.abs(D).sum(0)+np.abs(D).sum(1)>0)[0]
        u,v=loc(s.r0,c1),loc(s.r0,c2); U=np.array([u,v]); g=s.G[u,v]
        R=np.setdiff1d(idx,U)
        if len(R)==0: return abs(g)
        GUUi=np.array([[0,-1/g],[1/g,0]])
        H=s.G[np.ix_(R,R)]-s.G[np.ix_(R,U)]@GUUi@s.G[np.ix_(U,R)]
        sg,ld=np.linalg.slogdet(np.eye(len(R))+H@D[np.ix_(R,R)])
        return float(np.exp(0.5*ld)*abs(g))
if __name__=='__main__':
    import json,time
    b=json.load(open('../data/paper2_pair_bulk.json'))
    for t in ['0.3','0.5']:
        T0=time.time(); P=Pair(64,float(t)); ref=dict((int(r),c) for r,c in b[t])
        e=max(abs(P.C(32-(r+1)//2,32-(r+1)//2+r)-ref[r]) for r in range(1,33))
        print(t,'max abs err vs direct 64x64: %.1e  (%.1fs)'%(e,time.time()-T0))
