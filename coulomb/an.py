import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import collections
for N in [8,10,12,14,16]:
    W={}
    for l in open(f"W{N}.txt"):
        a,b,w=l.split(); W[(int(a),int(b))]=float(w)
    tot=sum(W.values()); C=N*N; n=C//2
    P=collections.defaultdict(float)
    for (a,b),w in W.items(): P[a]+=w/tot; P[b]+=w/tot
    rc=lambda x:divmod(x,N)
    adj=P[1]  # b1-like square
    # corner of unremoved corner diag neighbour (1,N-2)
    low=P[1*N+N-2]
    # screening: one defect in each half (r+c<N-1 vs >N-1)
    scr=0; sep=0; d_mean=0
    for (a,b),w in W.items():
        (ar,ac),(br,bc)=rc(a),rc(b)
        if (ar+ac<N-1 and br+bc>N-1) or (ar+ac>N-1 and br+bc<N-1): scr+=w/tot
        d_mean+=w/tot*((ar-br)**2+(ac-bc)**2)**.5
    # uniform comparison
    import itertools
    odd=[x for x in range(C) if sum(rc(x))%2]
    scr_u=sum(1 for a,b in itertools.combinations(odd,2) if (sum(rc(a))<N-1 and sum(rc(b))>N-1) or (sum(rc(a))>N-1 and sum(rc(b))<N-1))/(len(odd)*(len(odd)-1)/2)
    du=sum(((rc(a)[0]-rc(b)[0])**2+(rc(a)[1]-rc(b)[1])**2)**.5 for a,b in itertools.combinations(odd,2))/(len(odd)*(len(odd)-1)/2)
    # prob defect within Chebyshev radius R of a removed corner, R = N/4
    R=N//4
    near=sum(p for x,p in P.items() if max(rc(x))<R or min(rc(x))>N-1-R)/2
    near_u=sum(1 for x in odd if max(rc(x))<R or min(rc(x))>N-1-R)/len(odd)
    print(f"N={N:2d} unif={2/n*100:.3f}% adj={adj*100:.3f}% (x{adj/(2/n):.2f}) lowest-diag={low*100:.3f}% (x{low/(2/n):.2f}) "
          f"screen={scr:.3f} vs unif {scr_u:.3f}  <dist>/N={d_mean/N:.3f} vs {du/N:.3f}  near(N/4 box)={near:.3f} vs {near_u:.3f}")
    if N==16:
        for r in range(N-1,-1,-1):
            print(" ".join(f"{P[r*N+c]/(2/n):4.2f}" if (r+c)%2 else "  . " for c in range(N)))
