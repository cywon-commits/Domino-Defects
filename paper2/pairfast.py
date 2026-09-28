"""Fast pair correlation C(r) = Z(board-{u,v})/Z(board) for u=(r0,c1), v=(r0,c2) on an even board,
using a horizontal string between u and v (flip vertical edges (r0,c)-(r0+1,c), c1<c<c2, and
diagonals (r0,c)-(r0+1,c+1), c1<=c<c2).  K' = K + D with D supported on rows r0, r0+1.
Z(board-{u,v})/Z(board) = |Pf K'/Pf K| * |K'^{-1}(u,v)|, Pf K'/Pf K = ±sqrt(det(I + G_SS D_SS))."""
import numpy as np, scipy.linalg as sl, sys; sys.path.insert(0,'.')
from pftri import kmat
class Pair:
    def __init__(s, L, t, r0=None):
        s.L=L; s.t=t; s.r0 = L//2 if r0 is None else r0
        K=kmat(L,L,t); s.K=K
        S=np.r_[s.r0*L:(s.r0+2)*L]; s.S=S
        lu=sl.lu_factor(K, overwrite_a=False, check_finite=False)
        E=np.zeros((L*L,len(S))); E[S,np.arange(len(S))]=1
        s.G=sl.lu_solve(lu,E,check_finite=False)[S,:]   # G_SS
        s.KS=K[np.ix_(S,S)]; del lu
    def C(s, c1, c2):
        L=s.L; loc=lambda r,c:(r-s.r0)*L+c
        D=np.zeros_like(s.KS)
        for c in range(c1,c2+1):
            if c1<c<c2:
                i,j=loc(s.r0,c),loc(s.r0+1,c); D[i,j]-=2*s.KS[i,j]; D[j,i]-=2*s.KS[j,i]
            if c<c2 and c+1<L:
                i,j=loc(s.r0,c),loc(s.r0+1,c+1); D[i,j]-=2*s.KS[i,j]; D[j,i]-=2*s.KS[j,i]
        idx=np.where(np.abs(D).sum(0)+np.abs(D).sum(1)>0)[0]
        u,v=loc(s.r0,c1),loc(s.r0,c2); U=np.array([u,v])
        g=s.G[u,v]
        # K_{-U} (u,v deleted) is a valid Pfaffian orientation after the string; use
        # det(K_{-U}+D_R) = det K * g^2 * det(I + H_RR D_RR),  H = (K_{-U})^{-1} = G - G_.U G_UU^{-1} G_U.
        R=np.setdiff1d(idx,U)
        if len(R)==0: return abs(g)
        GUUi=np.array([[0,-1/g],[1/g,0]])          # inverse of [[0,g],[-g,0]]
        H=s.G[np.ix_(R,R)] - s.G[np.ix_(R,U)]@GUUi@s.G[np.ix_(U,R)]
        sign,ld=np.linalg.slogdet(np.eye(len(R))+H@D[np.ix_(R,R)])
        return float(np.exp(0.5*ld)*abs(g))

if __name__=='__main__':
    from pftri import logZ, logZ_multi
    for L,t in [(12,0.3),(12,0),(14,1.0)]:
        P=Pair(L,t); lz=logZ(L,L,t); err=0
        for c1 in range(0,L-1,3):
            for c2 in range(c1+1,L):
                ex=np.exp(logZ_multi(L,L,t,[(P.r0,c1),(P.r0,c2)])-lz); fa=P.C(c1,c2)
                err=max(err,abs(fa-ex)/max(ex,1e-300)) if ex>0 else max(err,fa)
        print(L,t,'max rel err %.1e'%err)
