"""Honeycomb flakes: edge probabilities from the Kasteleyn (biadjacency, all +1) inverse.
Used to test which boundaries give a FLAT liquid (all three orientations 1/3)."""
import numpy as np
s3=np.sqrt(3)
D=[np.array([0,1.0]),np.array([s3/2,-0.5]),np.array([-s3/2,-0.5])]   # A -> B bond vectors
def flake(R,kind):
    a1=np.array([s3,0]); a2=np.array([s3/2,1.5]); M=int(R/1.2)+4
    A=[];
    for m in range(-M,M+1):
        for n in range(-M,M+1): A.append(m*a1+n*a2)
    A=np.array(A); B=A+D[0]
    ang0 = 0 if kind=='armchair' else np.pi/6   # edge normals
    nrm=[np.array([np.cos(ang0+k*np.pi/3),np.sin(ang0+k*np.pi/3)]) for k in range(6)]
    c0=np.array([0,-1.0])   # hexagon (face) centre -> D6-symmetric, balanced flakes
    inside=lambda P: np.all([(P-c0)@v<=R for v in nrm],axis=0)
    A=A[inside(A)]; B=B[inside(B)]
    # prune dangling vertices repeatedly
    while True:
        key=lambda p:(round(p[0],4),round(p[1],4))
        Bset={key(b):i for i,b in enumerate(B)}; Aset={key(a):i for i,a in enumerate(A)}
        degA=np.array([sum(key(a+d) in Bset for d in D) for a in A]); degB=np.array([sum(key(b-d) in Aset for d in D) for b in B])
        if (degA>=2).all() and (degB>=2).all(): break
        A=A[degA>=2]; B=B[degB>=2]
    return A,B
def probs(A,B):
    key=lambda p:(round(p[0],4),round(p[1],4)); Bset={key(b):i for i,b in enumerate(B)}
    K=np.zeros((len(A),len(B))); E=[]
    for i,a in enumerate(A):
        for k,d in enumerate(D):
            j=Bset.get(key(a+d))
            if j is not None: K[i,j]=1; E.append((i,j,k))
    if len(A)!=len(B): return None,None
    Ki=np.linalg.inv(K)
    P=np.array([Ki[j,i] for i,j,k in E]); return E,P
if __name__=='__main__':
    for kind in ['armchair','zigzag']:
        for R in [10,18]:
            A,B=flake(R,kind); E,P=probs(A,B)
            if E is None: print(kind,R,'unbalanced',len(A),len(B)); continue
            k=np.array([e[2] for e in E]); pos=np.array([(A[e[0]]+B[e[1]])/2 for e in E])
            r=np.linalg.norm(pos-np.array([0,-1.0]),axis=1)/R
            out='%-8s R=%d sites=%d  '%(kind,R,2*len(A))
            for lo,hi in [(0,0.3),(0.3,0.6),(0.6,0.85),(0.85,1.01)]:
                m=(r>=lo)&(r<hi)
                out+=' r/R∈[%.2f,%.2f): '%(lo,hi)+'/'.join('%.3f'%P[m&(k==q)].mean() for q in range(3))
                out+=' frozen%%=%.0f'%(100*np.mean((P[m]<0.02)|(P[m]>0.98)))
            print(out)
