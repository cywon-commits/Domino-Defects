import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np, math
class Board:
    """Full N x N board, dense inverse Kasteleyn G(b,w). Defects removed via bordered determinant."""
    def __init__(s,N):
        s.N=N
        W=[(y,x) for y in range(N) for x in range(N) if (x+y)%2==1]; B=[(y,x) for y in range(N) for x in range(N) if (x+y)%2==0]
        s.wi={v:i for i,v in enumerate(W)}; s.bi={v:i for i,v in enumerate(B)}
        K=np.zeros((len(W),len(B)),complex)
        for (y,x),i in s.wi.items():
            for dy,dx,wt in ((0,1,1),(0,-1,1),(1,0,1j),(-1,0,1j)):
                j=s.bi.get((y+dy,x+dx))
                if j is not None: K[i,j]=wt
        s.G=np.linalg.inv(K)          # G[b,w]
    def string(s,v):
        y,x=v; N=s.N
        ks=range(y) if y<N/2 else range(y+1,N)
        xl=x if x<N-1 else x-1
        return [((k,xl),(k,xl+1)) for k in ks]
    def ratio(s,whites,blacks):
        """Z(board - whites - blacks)/Z(board). whites need strings; blacks must be on the boundary."""
        cnt={}
        for v in whites:
            for e in s.string(v): cnt[e]=cnt.get(e,0)^1
        edges=[e for e,c in cnt.items() if c]
        ew=[];eb=[]
        for a,b in edges:
            w,bl=(a,b) if (a[0]+a[1])%2==1 else (b,a); ew.append(s.wi[w]); eb.append(s.bi[bl])
        m=len(edges); k=len(whites)
        rows=eb+[s.bi[c] for c in blacks]; cols=ew+[s.wi[v] for v in whites]
        M=s.G[np.ix_(rows,cols)].copy()
        M[:m,:]*=-2; M[:m,:m]+=np.eye(m)
        return abs(np.linalg.det(M)) if m+k<300 else math.exp(np.linalg.slogdet(M)[1])
if __name__=="__main__":
    from kast import logcount
    N=24; Bd=Board(N); L0=logcount(N,[],[])
    for v1,v2 in [((5,6),(17,10)),((3,4),(20,21)),((12,13),(12,15)),((20,3),(21,6))]:
        r=Bd.ratio([v1,v2],[(0,0),(N-1,N-1)]); ref=math.exp(logcount(N,[(0,0),(N-1,N-1)],[v1,v2])-L0)
        print(v1,v2,r,ref)
