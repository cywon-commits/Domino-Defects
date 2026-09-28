"""Fit C(r) = Cinf + A r^-a e^{-ms r} + (-1)^r B r^-b e^{-ma r} on windows; compare ms, ma with 2kappa."""
import json, sys, numpy as np, warnings; from scipy.optimize import least_squares
sys.path.insert(0,'.'); from symbol_xi import kappa
warnings.filterwarnings('ignore')
d=json.load(open('../data/paper2_pair_big.json'))
def model(p,r):
    Ci,lA,a,ms,lB,b,ma,sA,sB=p[0],p[1],p[2],p[3],p[4],p[5],p[6],SA,SB
    return Ci + sA*np.exp(lA-a*np.log(r)-ms*r) + sB*(-1)**r*np.exp(lB-b*np.log(r)-ma*r)
def fit(r,C,k,fix_ab=None):
    global SA,SB
    Ci0=C[-1]; best=None
    for SA in [1,-1]:
        for SB in [1,-1]:
            p0=[Ci0,-2,1.0,2*k,-2,1.0,2*k]
            lo=[Ci0-1e-3,-40,-3,0.2*k,-40,-3,0.2*k]; hi=[Ci0+1e-3,5,8,8*k,5,8,8*k]
            def res(p):
                q=list(p)
                if fix_ab is not None: q[2]=fix_ab; q[5]=fix_ab
                return (model(q,r)-C)/(np.abs(C-Ci0)+1e-12)
            try: s=least_squares(res,p0,bounds=(lo,hi),xtol=1e-14,ftol=1e-14,max_nfev=20000)
            except Exception: continue
            if best is None or s.cost<best[0]: best=(s.cost,s.x.copy(),SA,SB)
    return best
if __name__=='__main__':
    for t in sys.argv[1:]:
        C=np.array(d[t]['C']); r=C[:,0]; C=C[:,1]; k=kappa(float(t),2001)
        tail=np.abs(C-C[-1])
        rmax=int(r[np.where(tail>1e-11*C[-1])[0].max()]) if (tail>1e-11).any() else 20
        print('t=%s  2kappa=%.4f  2artanh(t/sqrt(1-t^2))=%.4f   usable r<=%d'%(t,2*k,2*np.arctanh(float(t)/np.sqrt(1-float(t)**2)) if float(t)<0.707 else np.inf,rmax))
        for rmin in [4,8,12,16]:
            m=(r>=rmin)&(r<=rmax)
            if m.sum()<12: continue
            b=fit(r[m],C[m],k)
            if b is None: continue
            p=b[1]; print('   rmin=%2d  ms=%.4f (a=%.2f, A%s)   ma=%.4f (b=%.2f, B%s)   Cinf=%.13f  cost=%.1e'%(rmin,p[3],p[2],'+' if b[2]>0 else '-',p[6],p[5],'+' if b[3]>0 else '-',p[0],b[0]))
