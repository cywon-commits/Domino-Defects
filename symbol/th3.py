import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
exec(open('th2.py').read().split("for r in (1,3,5,7,9):")[0])
def bulk(r):
    def Gb(b,w): return G(abs(b[1]-w[1]),b[0]-w[0])
    edges=[]
    for k in range(r):
        A_,B_=(k,0),(k,1); ww,bl=(A_,B_) if (A_[0]+A_[1])%2 else (B_,A_); edges.append((ww,bl))
    v=(r,0); c0=(0,0); m=r; M=acb_mat(m+1,m+1)
    for j,(wj,bj) in enumerate(edges):
        for l,(wl,_) in enumerate(edges): M[j,l]=(1 if j==l else 0)-2*Gb(bj,wl)
        M[j,m]=-2*Gb(bj,v)
    for l,(wl,_) in enumerate(edges): M[m,l]=Gb(c0,wl)
    M[m,m]=Gb(c0,v); return abs(M.det())
for r in (1,3,5,7):
    print('bulk',r,bulk(r).str(15),(abs(D(r,lambda j,l:c(j-l-1)))/2).str(15))
