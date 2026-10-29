"""Infinite armchair strip |x| <= R (edges with normals 0 and 180 = classes + and -), periodic along y
with period 3. Exact bond probabilities by Fourier integration over k (trapezoid, periodic integrand).
A at m a1 + n a2, B = A + (0,1); bonds A -> A + d."""
import numpy as np
s3=np.sqrt(3); D=[np.array([0,1.0]),np.array([s3/2,-0.5]),np.array([-s3/2,-0.5])]
def build(R):
    M=int(R/s3)+4
    A=[];B=[]
    for n in range(2):                      # A rows y = 0, 1.5 in one period
        for m in range(-2*M,2*M+1):
            x=m*s3+n*s3/2
            if abs(x)<=R: A.append((x,1.5*n))
    for (x,y) in list(A):
        if True: B.append((x,y+1.0))
    A=np.array(A); B=np.array(B)
    key=lambda x,y:(round(x,4),round(y%3,4))
    while True:
        Bi={key(*b):j for j,b in enumerate(B)}; Ai={key(*a):i for i,a in enumerate(A)}
        dA=np.array([sum(key(*(a+d)) in Bi for d in D) for a in A]); dB=np.array([sum(key(*(b-d)) in Ai for d in D) for b in B])
        if (dA>=2).all() and (dB>=2).all(): break
        A=A[dA>=2]; B=B[dB>=2]
    Bi={key(*b):j for j,b in enumerate(B)}
    bonds=[]                                  # (i, j, s, dtype): A_i -> B_j shifted by s periods
    for i,a in enumerate(A):
        for t,d in enumerate(D):
            p=a+d; j=Bi.get(key(*p))
            if j is not None: bonds.append((i,j,int(np.floor((p[1]-B[j][1])/3+0.5)),t))
    return A,B,bonds
def probs(A,B,bonds,nk=256):
    nA=len(A); assert nA==len(B), (nA,len(B))
    ks=2*np.pi*(np.arange(nk)+0.5)/nk; P=np.zeros(len(bonds)); mindet=np.inf
    for k in ks:
        K=np.zeros((nA,nA),complex)
        for (i,j,s,t) in bonds: K[i,j]+=np.exp(1j*k*s)
        mindet=min(mindet,abs(np.linalg.det(K)) if nA<60 else np.min(np.abs(np.linalg.svd(K,compute_uv=False))))
        G=np.linalg.inv(K)
        for b,(i,j,s,t) in enumerate(bonds): P[b]+=(np.exp(1j*k*s)*G[j,i]).real
    return P/nk, mindet
if __name__=='__main__':
    import sys
    for R in [float(x) for x in sys.argv[1:]]:
        A,B,bonds=build(R)
        if len(A)!=len(B): print('R=%g unbalanced %d %d'%(R,len(A),len(B))); continue
        P,md=probs(A,B,bonds,128); P2,_=probs(A,B,bonds,256)
        vs=np.zeros(len(A))
        for (i,j,s,t),p in zip(bonds,P2): vs[i]+=p
        print('R=%g nA/period=%d  min singular value of K(k)=%.2e  vertex sums in [%.12f, %.12f]  nk 128 vs 256: %.1e'%(R,len(A),md,vs.min(),vs.max(),np.abs(P-P2).max()))
