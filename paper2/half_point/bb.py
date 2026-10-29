"""Basor-Bleher block Toeplitz representation K2(n) = 1/2 sqrt(det T_n(phi)) (arXiv:1610.08021, eqs 2.6-2.11)."""
import mpmath as mp
def coeffs(t, kmax, N=512, dps=60):
    mp.mp.dps=dps; t=mp.mpf(t)
    xs=[2*mp.pi*j/N for j in range(N)]
    vals=[]
    for x in xs:
        s,c=mp.sin(x),mp.cos(x)
        sig=1/((1-2*t*c+t*t)*mp.sqrt(t*t+s*s+s**4))
        p=lambda xx:(t*mp.cos(xx)+mp.sin(xx)**2)*(t-mp.e**(1j*xx))
        q=lambda xx:mp.sin(xx)*(1-2*t*mp.cos(xx)+t*t)
        vals.append([sig*p(x),sig*q(x),sig*q(-x),sig*p(-x)])
    co={}
    for k in range(-kmax,kmax+1):
        co[k]=[sum(vals[j][e]*mp.e**(-1j*k*xs[j]) for j in range(N))/N for e in range(4)]
    return co
def K2(co,n):
    M=mp.matrix(2*n,2*n)
    for a in range(n):
        for b in range(n):
            c=co[a-b]; M[2*a,2*b],M[2*a,2*b+1],M[2*a+1,2*b],M[2*a+1,2*b+1]=c
    return mp.sqrt(mp.re(mp.det(M)))/2
if __name__=='__main__':
    import sys
    t=sys.argv[1]; nmax=int(sys.argv[2]); dps=int(sys.argv[3])
    co=coeffs(t,nmax+1,dps=dps+10); mp.mp.dps=dps
    q=mp.sqrt(3)-2
    for n in range(1,nmax+1):
        k=K2(co,n); print(n, mp.nstr(k,25), '  closed forms (1-q^(n+s))/6 for s=0,1:', mp.nstr(k-(1-q**n)/6,3), mp.nstr(k-(1-q**(n+1))/6,3))
