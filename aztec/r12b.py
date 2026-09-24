import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
from kuo_lib import Z, aztec
import sympy as sp, pickle
x,y=sp.symbols('x y'); D={}
for n in range(2,9):
    A=aztec(n); ZA=2**(n*(n+1)//2)
    S1=[(n-k,k-1) for k in range(1,n+1)]; S2=[(n+l-1,l-1) for l in range(1,n+1)]
    T=[[Z(A-{w,b})*2**n//ZA for b in S2] for w in S1]
    assert all(Z(A-{w,b})*2**n%ZA==0 for w in S1 for b in S2)
    D[n]=T
    P=sum(T[k][l]*x**k*y**l for k in range(n) for l in range(n))
    print(n, sp.factor(P))
pickle.dump(D,open('r12tab.pkl','wb'))
