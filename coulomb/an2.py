import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import collections,itertools,math
for N in [8,12,16]:
    W={}
    for l in open(f"F{N}.txt"):
        a,b,w=l.split(); W[(int(a),int(b))]=float(w)
    tot=sum(W.values()); rc=lambda x:divmod(x,N)
    D=collections.defaultdict(float); U=collections.Counter()
    P=collections.defaultdict(float)
    for (a,b),w in W.items():
        (ar,ac),(br,bc)=rc(a),rc(b); d=abs(ar-br)+abs(ac-bc)
        D[d]+=w/tot; U[d]+=1; P[a]+=w/tot; P[b]+=w/tot
    nU=sum(U.values())
    adjP=D[1]; adjU=U[1]/nU
    md=sum(d*p for d,p in D.items()); mu=sum(d*c for d,c in U.items())/nU
    print(f"N={N}: P(adjacent)={adjP:.4f} vs uniform {adjU:.4f} (x{adjP/adjU:.1f}); P(d<=3)={sum(D[d] for d in (1,3)):.3f} vs {(U[1]+U[3])/nU:.3f}; mean dist/N={md/N:.3f} vs {mu/N:.3f}; median d={next(d for d in sorted(D) if sum(D[e] for e in D if e<=d)>=.5)}")
    if N==8:
        print(" d:", " ".join(f"{d}:{D[d]:.3f}/{U[d]/nU:.3f}" for d in sorted(D)))
        u=2/64
        for r in range(N-1,-1,-1): print(" ".join(f"{P[r*N+c]/u:4.2f}" for c in range(N)))
