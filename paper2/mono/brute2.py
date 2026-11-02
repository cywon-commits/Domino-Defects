import sys,numpy as np
Lx,Ly=int(sys.argv[1]),int(sys.argv[2]);J=float(sys.argv[3]);rem=set(int(a) for a in sys.argv[4:])
cells=[(x,y) for y in range(Ly) for x in range(Lx) if y*Lx+x not in rem]
def gen(free):
    if not free: yield (); return
    c=min(free,key=lambda p:(p[1],p[0])); x,y=c
    for q in [(x+1,y),(x,y+1)]:
        if q in free:
            for r in gen(free-{c,q}): yield ((c,q),)+r
def npl(dom):
    m={}
    for a,b in dom: m[a]=b;m[b]=a
    s=0
    for y in range(Ly-1):
        for x in range(Lx-1):
            a,b,c,d=(x,y),(x+1,y),(x,y+1),(x+1,y+1)
            if (m.get(a)==b and m.get(c)==d) or (m.get(a)==c and m.get(b)==d): s+=1
    return s
w=np.zeros(Lx*Ly); S=frozenset(cells)
for v in cells:
    w[v[1]*Lx+v[0]]=sum(np.exp(-J*npl(d)) for d in gen(S-{v}))
print(' '.join('%.8f'%x for x in w/w.sum()))
