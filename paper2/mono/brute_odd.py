import sys,numpy as np,functools
L=int(sys.argv[1]);J=float(sys.argv[2])
cells=[(r,c) for r in range(L) for c in range(L)]
def gen(free):
    if not free: yield (); return
    r,c=min(free)
    for q in [(r,c+1),(r+1,c)]:
        if q in free:
            for rest in gen(free-{(r,c),q}): yield (((r,c),q),)+rest
def npl(dom):
    m={}
    for a,b in dom: m[a]=b;m[b]=a
    s=0
    for r in range(L-1):
        for c in range(L-1):
            a,b,d,e=(r,c),(r,c+1),(r+1,c),(r+1,c+1)
            if (m.get(a)==b and m.get(d)==e) or (m.get(a)==d and m.get(b)==e): s+=1
    return s
w=np.zeros(L*L)
for v in cells:
    if (v[0]+v[1])%2: continue
    w[v[0]*L+v[1]]=sum(np.exp(-J*npl(d)) for d in gen(frozenset(cells)-{v}))
print(' '.join('%.6f'%x for x in w/w.sum()))
