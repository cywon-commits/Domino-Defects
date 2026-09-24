import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
exec(open('th2.py').read().split("for r in (1,3,5,7,9):")[0])
ctx.prec=500
cache={}
def cc(m):
    m=abs(m)
    if m not in cache: cache[m]=c(m)
    return cache[m]
for r in (1,3,5,7,9,15,21,31):
    Dr=D(r,lambda j,k: cc(j-k)-cc(j+k+2))
    mn=D(r,lambda j,k: cc(j-k-1)-cc(j+k+3))
    Dr1=D(r+1,lambda j,k: cc(j-k)-cc(j+k+2))
    print(r,'D_r',abs(Dr).str(10),'D_r+1',abs(Dr1).str(10),'minor/2',(abs(mn)/2).str(10),'p0',(abs(mn)/abs(Dr)).str(10))
