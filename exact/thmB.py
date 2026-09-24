import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np, math
from scipy.special import roots_legendre as leggauss
# ---------- quadrant inverse Kasteleyn via sine transform in x, closed form in y ----------
class Quadrant:
    def __init__(s, nq=4000):
        t,w=leggauss(nq)
        th=np.concatenate([(t+1)*np.pi/4, (t+1)*np.pi/4+np.pi/2]); wt=np.concatenate([w,w])*np.pi/4
        s.th,s.wt=th,wt
        c=2*np.cos(th); sg=np.sign(c); sq=np.sqrt(c*c+4)
        s.zin=1j*(c-sg*sq)/2; s.pref=sg/sq
    def F(s,nmax):   # F[n,theta] = g_theta(n), n=0..nmax
        n=np.arange(nmax+1)[:,None]
        return s.pref[None,:]*s.zin[None,:]**n
    def G(s, pts_b, pts_w):
        """matrix G_Q(b,w) for lists of (y,x) cells"""
        yb=np.array([p[0] for p in pts_b]); xb=np.array([p[1] for p in pts_b])
        yw=np.array([p[0] for p in pts_w]); xw=np.array([p[1] for p in pts_w])
        nmax=int(max(yb.max()+yw.max()+2, 1))
        def Tvec(P):
            T=np.empty(nmax+1,complex); blk=48; base=s.pref.copy()
            logz=np.log(s.zin)
            for a in range(0,nmax+1,blk):
                n=np.arange(a,min(a+blk,nmax+1))[:,None]
                T[a:a+len(n)]=(s.pref[None,:]*np.exp(n*logz[None,:]))@P
            return T
        xs=np.unique(np.concatenate([xb,xw]))
        S={x:np.sin(s.th*(x+1)) for x in xs}
        out=np.zeros((len(pts_b),len(pts_w)),complex)
        for x1 in np.unique(xb):
            for x2 in np.unique(xw):
                P=(2/np.pi)*s.wt*S[x1]*S[x2]
                T=Tvec(P)
                ib=np.where(xb==x1)[0]; iw=np.where(xw==x2)[0]
                D=np.abs(yb[ib][:,None]-yw[iw][None,:]); Sm=yb[ib][:,None]+yw[iw][None,:]+2
                out[np.ix_(ib,iw)]=T[D]-T[Sm]
        return out
def ratio(Gfun, v, corner=(0,0)):
    """Z(Q - corner - v)/Z(Q) for white interior defect v=(y,x); string = horizontal edges (k,x)-(k,x+1), k<y"""
    y,x=v
    edges=[]
    for k in range(y):
        a,b=(k,x),(k,x+1)
        w,bl=(a,b) if (a[0]+a[1])%2==1 else (b,a)
        edges.append((w,bl))
    W=[e[0] for e in edges]; B=[e[1] for e in edges]
    m=len(edges)
    GB=Gfun(B+[corner], W+[v])          # rows: blacks (path + corner), cols: whites (path + v)
    Mx=np.zeros((m+1,m+1),complex)
    Mx[:m,:m]=np.eye(m)-2*GB[:m,:m]      # I + A,  A_jl = -2 K(w_j,b_j) G(b_j,w_l), K=1 on horizontal edges
    Mx[:m,m]=-2*GB[:m,m]                 # q
    Mx[m,:m]=GB[m,:m]                    # s
    Mx[m,m]=GB[m,m]                      # t
    sign,ld=np.linalg.slogdet(Mx)
    return math.exp(ld)
