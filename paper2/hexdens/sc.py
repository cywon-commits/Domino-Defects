"""Conformal map of the regular hexagon (corners at 30+60k deg, apothem a) and the Coulomb-gas density."""
import numpy as np
from scipy.special import beta, poch, gammaln
NT=8000
n=np.arange(NT); an=np.exp(gammaln(n+1/3)-gammaln(1/3)-gammaln(n+1))     # (1/3)_n/n!
V=beta(1/6,2/3)/6                                                         # int_0^1 (1-x^6)^(-1/3)
ROT=np.exp(1j*np.pi/6)
def F(z):   # sum a_n z^(6n+1)/(6n+1)  (vectorised over z, Horner in w=z^6)
    w=z**6; coef=an/(6*n+1); s=np.zeros_like(z)
    for c in coef[::-1]: s=s*w+c
    return z*s
def F_horner_fast(z):
    w=z**6; coef=(an/(6*n+1))
    # chunked power sum for speed
    out=np.zeros_like(z); p=np.ones_like(z)
    for k in range(0,NT,2000):
        c=coef[k:k+2000]; P=p[:,None]*w[:,None]**np.arange(len(c))[None,:]
        out+=P@c; p=p*w**len(c)
    return z*out
def inv(Z,a):
    """zeta in unit disk with hexagon point Z (apothem a)."""
    Rc=a/np.cos(np.pi/6); C=Rc/V; target=Z/(C*ROT)       # target in the hexagon with corners at 0,60,..., circumradius V
    zeta=target/V*0.95
    for it in range(60):
        Fz=F_horner_fast(zeta); dF=(1-zeta**6)**(-1/3)
        step=(Fz-target)/dF; zeta=zeta-step
        r=np.abs(zeta); big=r>0.999999; zeta[big]=zeta[big]/r[big]*0.999999
        if np.max(np.abs(step))<1e-13: break
    return zeta, C
def density(Z,a,q,e=1,g=0.5,pre=None):
    """q: charges at corners 30,90,150,210,270,330 deg. Returns log rho (up to a constant)."""
    zeta,C=inv(Z,a) if pre is None else pre
    al=np.pi/6; u=np.exp(-1j*al)
    f=1j*(1+u*zeta)/(1-u*zeta)
    dfdzeta=2j*u/(1-u*zeta)**2; dzdzeta=C*ROT*(1-zeta**6)**(-1/3)
    fp=dfdzeta/dzdzeta
    zc=np.exp(1j*np.pi*np.arange(6)/3)              # disk prevertices of corners at 30+60k (after ROT)
    fc=1j*(1+u*zc)/(1-u*zc)
    lr=(g/2)*np.log(np.abs(fp))+(g/2)*np.log(np.imag(f))
    for k in range(6): lr+=2*g*e*q[k]*np.log(np.abs(f-fc[k]))
    return lr, zeta
if __name__=='__main__':
    # self-test: corners map to |zeta|=1 at the right angles, and the boundary midpoint
    a=10.0; Rc=a/np.cos(np.pi/6)
    pts=np.array([0.95*Rc*np.exp(1j*np.radians(30+60*k)) for k in range(6)]+[0.999*a+0j, 0.5*a*1j])
    z,_=inv(pts,a); print(np.round(np.abs(z),4)); print(np.round(np.degrees(np.angle(z)),2))
    zt=F_horner_fast(z)*(Rc/V)*ROT; print('roundtrip max err',np.abs(zt-pts).max())
