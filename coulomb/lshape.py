import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import math, subprocess
from kreg import logcount_region
def L(N): n=N//2; return {(r,c) for r in range(N) for c in range(N) if not (r>=n and c>=n)}
# verify vs DP on N=8
N=8; R=L(N); out=[(r,c) for r in range(N) for c in range(N) if (r,c) not in R]
for rem,dfc in [((0,0),(4,3)),((0,0),(1,2)),((0,7),(3,3)),((0,7),(2,2))]:
    args=[str(N),str(N)]+[str(x) for v in out+[rem,dfc] for x in v]
    dp=float(subprocess.check_output(["./cntf"]+args)); lc=logcount_region(R-{rem},[dfc])
    print("check",rem,dfc,dp,round(math.exp(lc)))
print("N  attractive(white defect, white concave corner)  repulsive(black defect)")
prev=None
for N in (32,64,128,256):
    n=N//2; R=L(N)
    a=logcount_region(R-{(0,0)},[(n,n-1)])-logcount_region(R-{(0,0)},[(N//4,N//4+1)])
    b=logcount_region(R-{(0,N-1)},[(n-1,n-1)])-logcount_region(R-{(0,N-1)},[(N//4,N//4)])
    s=""
    if prev: s=f"  slopes: {-(a-prev[0])/math.log(2):+.3f}  {-(b-prev[1])/math.log(2):+.3f}"
    print(N, f"{math.exp(a):.4f}", f"{math.exp(b):.4f}", s, flush=True); prev=(a,b)
