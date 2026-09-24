import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
# Exact plane inverse-Kasteleyn G(x,y) = A + B/pi  (A,B Gaussian rationals), weights: horizontal 1, vertical i.
from flint import fmpq, arb, acb, acb_mat, ctx
class V:  # a + b/pi, a,b complex rationals (re,im)
    __slots__=("ar","ai","br","bi")
    def __init__(s,ar=0,ai=0,br=0,bi=0): s.ar,s.ai,s.br,s.bi=fmpq(ar),fmpq(ai),fmpq(br),fmpq(bi)
    def __add__(s,o): return V(s.ar+o.ar,s.ai+o.ai,s.br+o.br,s.bi+o.bi)
    def __sub__(s,o): return V(s.ar-o.ar,s.ai-o.ai,s.br-o.br,s.bi-o.bi)
    def neg(s): return V(-s.ar,-s.ai,-s.br,-s.bi)
    def muli(s): return V(-s.ai,s.ar,-s.bi,s.br)       # multiply by i
    def scale(s,q): return V(s.ar*q,s.ai*q,s.br*q,s.bi*q)
    def toacb(s,pi): return acb(arb(s.ar)+arb(s.br)/pi, arb(s.ai)+arb(s.bi)/pi)
ZERO=V()
def plane_table(Xmax,Ymax):
    """G[x][y] for 0<=x<=Xmax, 0<=y<=Ymax (use symmetry for negatives)."""
    Ytot=Ymax+Xmax+2
    # J_n = p_n*pi + q_n for odd n
    p={1:fmpq(1,4)}; q={1:fmpq(0)}; p[-1]=fmpq(0); q[-1]=fmpq(0)
    for n in range(1,Ytot+2,2):
        # J_{n+2} = (6n J_n - (n-1) J_{n-2} - 4)/(n+1)
        pm=p.get(n-2,fmpq(0)); qm=q.get(n-2,fmpq(0))
        p[n+2]=(6*n*p[n]-(n-1)*pm)/(n+1); q[n+2]=(6*n*q[n]-(n-1)*qm-4)/(n+1)
    col0=[]
    for n in range(Ytot+2):
        if n%2==0: col0.append(ZERO); continue
        # (-i)^n (p + q/pi)
        ph=n%4  # (-i)^1=-i, (-i)^3=i
        v=V(0,-p[n],0,-q[n]) if ph==1 else V(0,p[n],0,q[n])
        col0.append(v)
    g=lambda col,y: col[abs(y)]
    col1=[]
    for y in range(Ytot+1):
        s=(g(col0,y+1)+g(col0,y-1)).muli()
        d=V(1) if y==0 else ZERO
        col1.append((d-s).scale(fmpq(1,2)))
    cols=[col0,col1]
    for x in range(1,Xmax):
        prev,cur=cols[x-1],cols[x]; L=len(cur)-1
        new=[]
        for y in range(L):
            s=(g(cur,y+1)+g(cur,y-1)).muli()
            new.append((prev[y]+s).neg())
        cols.append(new)
    return cols
