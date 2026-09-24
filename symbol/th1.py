import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
from exactG import *
from flint import arb, acb, acb_mat, ctx
ctx.prec=400
R=30; cols=plane_table(3,4*R+10); pi=arb.pi()
G=lambda dx,dy: cols[abs(dx)][abs(dy)].toacb(pi)
def a(m): return G(1,m) if m%2==0 else G(0,m)
print('a(0..6)',[a(m).str(8) for m in range(7)])
def halfratio(r):
    # exact_half construction
    def Gh(b,w):
        dx=abs(b[1]-w[1]); return G(dx,b[0]-w[0])-G(dx,b[0]+w[0]+2)
    edges=[]
    for k in range(r):
        A_,B_=(k,0),(k,1); ww,bl=(A_,B_) if (A_[0]+A_[1])%2 else (B_,A_); edges.append((ww,bl))
    v=(r,0); c=(0,0); m=r
    M=acb_mat(m+1,m+1)
    for j,(wj,bj) in enumerate(edges):
        for l,(wl,_) in enumerate(edges):
            M[j,l]=(1 if j==l else 0)-2*Gh(bj,wl)
        M[j,m]=-2*Gh(bj,v)
    for l,(wl,_) in enumerate(edges): M[m,l]=Gh(c,wl)
    M[m,m]=Gh(c,v)
    return abs(M.det())
def TH(r,s,eps,shift=0):
    M=acb_mat(r,r)
    for j in range(r):
        for l in range(r):
            M[j,l]=(1 if j==l else 0)-2*a(j-l+shift)+eps*2*a(j+l+s)
    return M.det()
for r in (1,3,5,7):
    h=halfratio(r)
    print('r',r,'half',h.str(15))
    for s in range(0,4):
        for eps in (1,-1):
            for rr in (r-1,r,r+1):
                if rr<1: continue
                d=TH(rr,s,eps)
                print('    s',s,'eps',eps,'size',rr, abs(d).str(12))
