import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import itertools, functools
N=8
def count(holes):
    # broken-profile DP, row-major cells
    blocked=set(holes)
    @functools.lru_cache(None)
    def go(pos,mask):
        if pos==N*N: return 1 if mask==0 else 0
        r,c=divmod(pos,N)
        if (mask&1) or (r,c) in blocked:
            if (mask&1) and (r,c) in blocked: return 0
            return go(pos+1,mask>>1)
        t=0
        # vertical: occupy (r+1,c)
        if r+1<N and (r+1,c) not in blocked: t+=go(pos+1,(mask>>1)|(1<<(N-1)))
        # horizontal
        if c+1<N and (r,c+1) not in blocked and not (mask&2): t+=go(pos+1,(mask>>1)|1)
        return t
    return go(0,0)
corners=[(0,0),(7,7)]
print("full board",count([]))
free=[(r,c) for r in range(N) for c in range(N) if (r+c)%2==1]  # colour opposite to removed corners
tot=0; W={}
for a,b in itertools.combinations(free,2):
    w=count(corners+[a,b]); W[(a,b)]=w; tot+=w
print("total configs",tot)
P={s:0 for s in free}
for (a,b),w in W.items(): P[a]+=w/tot; P[b]+=w/tot
for r in range(N-1,-1,-1):
    print(" ".join(f"{P[(r,c)]*100:5.2f}" if (r,c) in P else "  -- " for c in range(N)))
print("uniform per square:",2/32*100)
import collections
D=collections.defaultdict(float);C=collections.Counter()
for (a,b),w in W.items():
    d=abs(a[0]-b[0])+abs(a[1]-b[1]); D[d]+=w/tot; C[d]+=1
for d in sorted(D): print(d, round(D[d],4), round(C[d]/len(W),4))
best=sorted(W.items(),key=lambda x:-x[1])[:5]; worst=sorted(W.items(),key=lambda x:x[1])[:5]
print(best);print(worst)
