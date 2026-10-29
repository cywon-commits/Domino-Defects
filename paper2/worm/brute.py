import sys, numpy as np, itertools
W,H,J,xt,yt=int(sys.argv[1]),int(sys.argv[2]),float(sys.argv[3]),int(sys.argv[4]),int(sys.argv[5])
idx=lambda x,y: y*W+(x%W)
N=W*H; tail=idx(xt,yt)
def nbrs(i):
    x,y=i%W,i//W; o=[idx(x+1,y),idx(x-1,y)]
    if y+1<H:o.append(idx(x,y+1))
    if y>0:o.append(idx(x,y-1))
    return o
def matchings(free):
    if not free: yield {}; return
    i=min(free)
    for j in nbrs(i):
        if j in free and j!=i:
            for m in matchings(free-{i,j}):
                mm=dict(m); mm[i]=j; mm[j]=i; yield mm
def nplaq(m):
    s=0
    for y in range(H-1):
        for x in range(W):
            a,b,c,d=idx(x,y),idx(x+1,y),idx(x,y+1),idx(x+1,y+1)
            if (m.get(a)==b and m.get(c)==d) or (m.get(a)==c and m.get(b)==d): s+=1
    return s
hist=np.zeros(N)
for h in range(N):
    if h==tail or ((h%W)+(h//W))%2==((xt+yt)%2): continue
    free=frozenset(range(N))-{tail,h}
    seen=set(); 
    for m in matchings(free):
        key=tuple(sorted((a,b) for a,b in m.items() if a<b))
        if key in seen: continue
        seen.add(key); hist[h]+=np.exp(-J*nplaq(m))
print(' '.join('%.6f'%v for v in hist/hist.sum()))
