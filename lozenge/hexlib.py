import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import flint, math
def hexagon(a,b,c):
    """point hexagon: 0<=i<=a+b? use parallelogram cut. sides: bottom (j=0) length a, right (i=A) length b, ... """
    # points with 0<=j<=b+c, 0<=i<=a+c, c<=i+j<=a+b+c  -> hexagon sides: bottom a, right(i=a+c) b, top-right c, top a, left b, bottom-left c
    A=a+c; B=b+c
    inside=lambda i,j: 0<=i<=A and 0<=j<=B and c<=i+j<=a+b+c
    U=[];D=[]
    for i in range(-1,A+1):
        for j in range(-1,B+1):
            if inside(i,j) and inside(i+1,j) and inside(i,j+1): U.append((i,j))
            if inside(i+1,j) and inside(i,j+1) and inside(i+1,j+1): D.append((i,j))
    return set(U),set(D)
def Kmat(U,D):
    U=sorted(U); D=sorted(D); iu={u:k for k,u in enumerate(U)}; idd={d:k for k,d in enumerate(D)}
    M=flint.fmpz_mat(len(U),len(D))
    for (i,j),r in iu.items():
        for d in ((i,j),(i-1,j),(i,j-1)):
            if d in idd: M[r,idd[d]]=1
    return M,iu,idd
def Z(U,D):
    if len(U)!=len(D): return 0
    M,_,_=Kmat(U,D); return abs(int(M.det()))
def macmahon(a,b,c):
    from fractions import Fraction
    p=Fraction(1)
    for i in range(1,a+1):
        for j in range(1,b+1):
            for k in range(1,c+1): p*=Fraction(i+j+k-1,i+j+k-2)
    return p
