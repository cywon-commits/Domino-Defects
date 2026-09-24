import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import math
for N in [4,6,8,10,12,14,16]:
    rows=[l.split() for l in open(f"O{N}.txt")]
    W={int(b):(int(hi)<<64)+int(lo) for b,_,hi,lo in rows}
    tot=sum(W.values()); n=len(W); u=1/n
    P={b:w/tot for b,w in W.items()}
    rc=lambda x:divmod(x,N)
    adj=P[1]; far=P[(N-1)*N+N-2]  # next to opposite corner h8 (h7 / g8)
    oc1=P[N-1] if (N-1)%2 else None  # h1 corner is white? (0,N-1): parity N-1 odd -> white
    md=sum(p*math.hypot(*rc(b)) for b,p in P.items())/N; mu=sum(math.hypot(*rc(b)) for b in P)/n/N
    print(f"N={N:2d} total={tot} adj={adj/u:.2f}x  h1/a8 corner={P[N-1]/u:.2f}x  near-opposite-corner={far/u:.2f}x  min={min(P.values())/u:.2f}x max={max(P.values())/u:.2f}x  <dist from a1>/N={md:.3f} (unif {mu:.3f})")
    if N==8:
        # integer weights
        for r in range(N-1,-1,-1): print("   "+" ".join(f"{W[r*N+c]:>9d}" if (r*N+c) in W else "        ." for c in range(N)))
        for r in range(N-1,-1,-1): print("   "+" ".join(f"{P[r*N+c]/u:5.2f}" if (r*N+c) in W else "    ." for c in range(N)))
    if N==16:
        for r in range(N-1,-1,-1): print("   "+" ".join(f"{P[r*N+c]/u:4.2f}" if (r*N+c) in W else "   ." for c in range(N)))
