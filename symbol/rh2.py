import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
from exactG import *
from flint import arb, acb, acb_mat, ctx
import sys, time
Rmax=int(sys.argv[1]); ctx.prec=int(Rmax*3.4)+200
cols=plane_table(1,2*Rmax+10); pi=arb.pi()
G=lambda dx,dy: cols[abs(dx)][abs(dy)].toacb(pi)
c=[ (1 if m==0 else 0)-2*(G(1,m) if m%2==0 else G(0,m)) for m in range(2*Rmax+8)]
cc=lambda m: c[abs(m)]
def det(n,f):
    M=acb_mat(n,n)
    for j in range(n):
        for k in range(n): M[j,k]=f(j,k)
    return M.det()
for r in [int(x) for x in sys.argv[2:]]:
    t=time.time()
    Dr=abs(det(r,lambda j,k: cc(j-k)-cc(j+k+2)))
    mn=abs(det(r,lambda j,k: cc(j-k-1)-cc(j+k+3)))
    print(r, Dr.str(25,radius=False), (mn/Dr).str(25,radius=False), (mn/2).str(25,radius=False), '%.0fs'%(time.time()-t), flush=True)
