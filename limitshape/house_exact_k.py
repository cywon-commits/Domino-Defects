import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import sys, math, json, time
from flint import fmpz_mat, fmpq
def house(N,k1,k2): return {(y,x) for y in range(N) for x in range(N) if (x+y)<=2*N-2-k1 and (y-x)<=N-1-k2}
out={}
FR=float(sys.argv[1])
for N in [int(a) for a in sys.argv[2:]]:
    t=time.time(); k=int(round(FR*N)); k-=k%2; cells=house(N,k,k)
    W=sorted(v for v in cells if sum(v)%2); B=sorted(v for v in cells if sum(v)%2==0)
    wi={v:i for i,v in enumerate(W)}; bi={v:i for i,v in enumerate(B)}
    M=fmpz_mat(len(W),len(B))
    for (y,x),i in wi.items():
        for dy,dx in ((0,1),(0,-1),(1,0),(-1,0)):
            j=bi.get((y+dy,x+dx))
            if j is not None: M[i,j]= 1 if dy==0 else (1 if x%2==0 else -1)
    e=fmpz_mat(len(B),1); e[bi[(0,0)],0]=1
    y=M.transpose().solve(e)
    bnd=lambda v: any((v[0]+a,v[1]+b) not in cells for a,b in ((1,0),(-1,0),(0,1),(0,-1)))
    prof={}
    for v in W:
        if bnd(v):
            q=abs(y[wi[v],0]); prof[str(v)]=(math.log(int(q.p))-math.log(int(q.q))) if q!=0 else None
    out[N]=prof
    stair=sorted([v for v in W if bnd(v) and v[0]-v[1]>=N-2-k], key=lambda v:v[1])
    vals=[(v[1]/N, prof[str(v)]/N) for v in stair if prof[str(v)] is not None]
    best=max(vals,key=lambda p:p[1]); pm=max(p for p in prof.values() if p is not None); btot=(pm+math.log(sum(math.exp(p-pm) for p in prof.values() if p is not None)))/N
    print(f"N={N}: boundary-weight total log/N={btot:.4f}; staircase max (1/N)ln w={best[1]:.4f} at x/N={best[0]:.3f}  ({time.time()-t:.0f}s)",flush=True)
json.dump({str(k):v for k,v in out.items()},open("house_exact_k%s.json"%sys.argv[1],"w"))
