import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
exec(open('th1.py').read().split("for r in (1,3,5,7):")[0])
def c(m): return (1 if m==0 else 0)-2*a(m)
def D(r,f):
    M=acb_mat(r,r)
    for j in range(r):
        for l in range(r): M[j,l]=f(j,l)
    return M.det()
for r in (1,3,5,7,9):
    h=halfratio(r)
    d=D(r,lambda j,l: c(j-l-1)-c(j+l+3))
    print(r,h.str(15),(abs(d)/2).str(15))
