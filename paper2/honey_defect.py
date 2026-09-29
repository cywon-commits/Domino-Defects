"""Task 2: forced defect on A-centred armchair honeycomb flakes (n_A = n_B + 1 -> one A defect).
Z(flake - v) for an A site v via sparse log|det| of the biadjacency matrix with a string:
flip every bond crossed by the downward half-line from v (offset eps to avoid vertices)."""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
s3=np.sqrt(3); D=[np.array([0,1.0]),np.array([s3/2,-0.5]),np.array([-s3/2,-0.5])]
def flakeA(R, c0=np.zeros(2)):
    a1=np.array([s3,0]); a2=np.array([s3/2,1.5]); M=int(R/1.2)+4
    A=np.array([m*a1+n*a2 for m in range(-M,M+1) for n in range(-M,M+1)]); B=A+D[0]
    nrm=[np.array([np.cos(k*np.pi/3),np.sin(k*np.pi/3)]) for k in range(6)]   # armchair edges
    ins=lambda P: np.all([(P-c0)@v<=R for v in nrm],axis=0)
    A=A[ins(A)]; B=B[ins(B)]
    key=lambda p:(round(p[0],4),round(p[1],4))
    while True:
        Bs={key(b) for b in B}; As={key(a) for a in A}
        dA=np.array([sum(key(a+d) in Bs for d in D) for a in A]); dB=np.array([sum(key(b-d) in As for d in D) for b in B])
        if (dA>=2).all() and (dB>=2).all(): break
        A=A[dA>=2]; B=B[dB>=2]
    return A,B
def bonds(A,B):
    key=lambda p:(round(p[0],4),round(p[1],4)); Bi={key(b):j for j,b in enumerate(B)}
    return [(i,Bi[key(a+d)]) for i,a in enumerate(A) for d in D if key(a+d) in Bi]
def crosses(p,q,x0,y0):
    # does segment p-q cross the half-line x=x0, y<y0 ?
    if (p[0]-x0)*(q[0]-x0)>=0: return False
    y=p[1]+(q[1]-p[1])*(x0-p[0])/(q[0]-p[0]); return y<y0
class Flake:
    def __init__(s,R):
        s.A,s.B=flakeA(R); s.E=bonds(s.A,s.B); s.R=R
    def logZ_minus(s,v):
        A,B=s.A,s.B; x0=A[v][0]+1e-3*np.pi; y0=A[v][1]
        I=[];J=[];V=[]
        for i,j in s.E:
            if i==v: continue
            w=-1.0 if crosses(A[i],B[j],x0,y0) else 1.0
            I.append(i if i<v else i-1); J.append(j); V.append(w)
        M=sp.csc_matrix((V,(I,J)),shape=(len(A)-1,len(B)))
        lu=spl.splu(M,permc_spec='COLAMD')
        return np.sum(np.log(np.abs(lu.U.diagonal())))
def brute(A,B,E,v):
    import functools
    adj={}
    for i,j in E:
        if i!=v: adj.setdefault(j,[]).append(i)
    Bs=tuple(range(len(B)))
    @functools.lru_cache(None)
    def f(usedA,k):
        if k==len(B): return 1
        return sum(f(usedA|(1<<i),k+1) for i in adj.get(k,[]) if not usedA>>i&1)
    return f(0,0)
if __name__=='__main__':
    for R in [3.5,4.5,5.2]:
        F=Flake(R); nA,nB=len(F.A),len(F.B)
        if nA!=nB+1: print(R,'nA-nB=',nA-nB); continue
        bad=0
        for v in range(nA):
            z=brute(F.A,F.B,F.E,v); p=np.exp(F.logZ_minus(v)) if z>0 else 0
            if z>0 and abs(p-z)>1e-6*z: bad+=1
        print('R=%.1f nA=%d  all %d defect sites checked, mismatches=%d'%(R,nA,nA,bad),flush=True)
