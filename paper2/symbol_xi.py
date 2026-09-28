"""Exact fermion decay rate along x from the Kasteleyn symbol (2-column unit cell).
K_AB(z,w) = 1 - 1/z - t w - t/(z w),  det K = -(w-1/w)^2 + K_AB(z,w) K_AB(1/z,1/w).
For real phi (w=e^{i phi}) solve det=0 for z; rate per cell = min |ln|z||; per site = /2."""
import numpy as np, sys; sys.path.insert(0,'.')
def kappa(t, nphi=4001):
    best=np.inf
    for phi in np.linspace(-np.pi,np.pi,nphi):
        w=np.exp(1j*phi)
        # K_AB(z,w)=a + b/z, a=1-t w, b=-(1+t/w); K_AB(1/z,1/w)=a' + b' z, a'=1-t/w, b'=-(1+t w)
        a,b=1-t*w,-(1+t/w); ap,bp=1-t/w,-(1+t*w); c=-(w-1/w)**2
        # (a+b/z)(ap+bp z)+c=0 -> a bp z^2 + (a ap + b bp + c) z + b ap = 0
        roots=np.roots([a*bp, a*ap+b*bp+c, b*ap])
        best=min(best,min(abs(np.log(abs(r))) for r in roots))
    return best/2   # per lattice site
if __name__=='__main__':
    # check against direct inverse of K on a finite board: G(u,v) along the middle row
    from pftri import kmat
    L=80; t=0.3; K=kmat(L,L,t); G=np.linalg.inv(K)
    r0=L//2; u=r0*L+L//4
    g=np.array([abs(G[u,u+d]) for d in range(1,41)])
    print('t=0.3 exact kappa/site = %.5f'%kappa(t),' artanh(t)=%.5f'%np.arctanh(t))
    print('finite-board |G| local rate (odd d):',' '.join('%.4f'%(np.log(g[d-1]/g[d+1])/2) for d in range(1,36,4)))
    for t in [0.05,0.1,0.2,0.3,0.5,0.7,1.0]:
        k=kappa(t); print('t=%.2f  kappa=%.5f  artanh(t)=%.5f  xi_f=1/kappa=%.3f  2kappa=%.4f'%(t,k,np.arctanh(t),1/k,2*k))
