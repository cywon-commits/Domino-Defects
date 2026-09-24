import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
from kuo_lib import Z, aztec, math
import sys
def setup(n):
    A=aztec(n); b=[(n-j,j) for j in range(n+1)]
    SW={v for v in A if v[0]+v[1]<=n}; NE={v for v in A if v[0]+v[1]>=3*n-2}
    return A,b,SW,NE,A-SW-NE
nb=lambda v:[(v[0]+d,v[1]+e) for d,e in((0,1),(0,-1),(1,0),(-1,0))]
def strip_forced(R):
    R=set(R); tiles=[]
    while True:
        ch=False
        for v in list(R):
            if v not in R: continue
            N=[u for u in nb(v) if u in R]
            if len(N)==1: R-={v,N[0]}; tiles.append((v,N[0])); ch=True
            elif len(N)==0: return None,tiles
        if not ch: return R,tiles
for n in range(3,9):
    A,b,SW,NE,I=setup(n)
    c=(1,n); d=(n,1)
    assert c in I and d in I
    ok1=all(Z((I|{b[m]})-{c})==2**((n-1)*(n-2)//2)*math.comb(n-1,m) for m in range(n+1))
    ok2=all(Z((I|{b[m]})-{d})==2**((n-1)*(n-2)//2)*math.comb(n-1,m-1) if m>0 else Z((I|{b[0]})-{d})==0 for m in range(n+1))
    zcd=Z(I-{c,d})
    # forced reduction of I - c
    R,t=strip_forced(I-{c})
    # compare with AD_{n-1} minus SW layer, translated
    A1,b1,SW1,NE1,I1=setup(n-1)
    T=A1-SW1
    shifts=[(dy,dx) for dy in range(-2,3) for dx in range(-2,3) if {(y+dy,x+dx) for y,x in T}==R]
    print(n, ok1, ok2, zcd==2**((n-1)*(n-2)//2), len(t), shifts)
print('cd check')
for n in range(3,11):
    A,b,SW,NE,I=setup(n); c=(1,n); d=(n,1)
    R,t=strip_forced(I-{c,d}); A2=aztec(n-2)
    sh=[(dy,dx) for dy in range(-3,4) for dx in range(-3,4) if {(y+dy,x+dx) for y,x in A2}==R]
    Rd,td=strip_forced(I-{d}); A1,b1,SW1,_,_=setup(n-1)
    shd=[(dy,dx) for dy in range(-2,3) for dx in range(-2,3) if {(y+dy,x+dx) for y,x in A1-SW1}==Rd]
    print(n,len(t),sh,len(td),shd)
    # final R11b check via Kuo formula
    ok=all(Z(I|{b[s],b[s2]})==2**((n-1)*(n-2)//2)*(math.comb(n-1,s2-1)*math.comb(n-1,s)-math.comb(n-1,s2)*math.comb(n-1,s-1) if s>0 else math.comb(n-1,s2-1)) for s in range(n+1) for s2 in range(s+1,n+1)) if n<=8 else None
    print('  R_{s,s2} formula',ok)
