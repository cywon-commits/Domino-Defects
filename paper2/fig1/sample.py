"""Uniform sample of dimer coverings of the 8x8 board minus corner (0,0) with one forced white monomer."""
import numpy as np, networkx as nx, random, json
L=8; rng=random.Random(11)
cells=[(r,c) for r in range(L) for c in range(L) if (r,c)!=(0,0)]
def nb(p): r,c=p; return [q for q in [(r+1,c),(r-1,c),(r,c+1),(r,c-1)] if q in S]
S=set(cells)
# initial: remove a white cell v, perfect matching of the rest
v=(3,4)
G=nx.Graph(); G.add_nodes_from(S-{v})
for p in S-{v}:
    for q in nb(p):
        if q!=v: G.add_edge(p,q)
M=nx.max_weight_matching(G,maxcardinality=True); mt={}
for a,b in M: mt[a]=b; mt[b]=a
assert len(mt)==len(S)-1
mono=v
for step in range(200000):
    if rng.random()<0.5:   # monomer move
        o=nb(mono); j=rng.choice(o); k=mt[j]; zk=len(nb(k))
        if rng.random()<len(o)/zk:
            mt[mono]=j; mt[j]=mono; del mt[k]; mono=k
    else:                  # plaquette rotation
        r,c=rng.randrange(L-1),rng.randrange(L-1)
        a,b,d,e=(r,c),(r,c+1),(r+1,c),(r+1,c+1)
        if not all(x in mt for x in (a,b,d,e)): continue
        if mt[a]==b and mt[d]==e: mt[a]=d;mt[d]=a;mt[b]=e;mt[e]=b
        elif mt[a]==d and mt[b]==e: mt[a]=b;mt[b]=a;mt[d]=e;mt[e]=d
    if step>150000 and 2<=mono[0]<=5 and 2<=mono[1]<=6: break
dom=sorted({tuple(sorted((a,b))) for a,b in mt.items()})
json.dump({'dominoes':dom,'mono':mono},open('sample.json','w')); print('monomer at',mono,'dominoes',len(dom))
