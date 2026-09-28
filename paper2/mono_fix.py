import json,sys,numpy as np,warnings; from scipy.optimize import least_squares
sys.path.insert(0,'.'); from symbol_xi import kappa; warnings.filterwarnings('ignore')
d=json.load(open('../data/paper2_pair_big.json'))
def run(t,rmin,ms_fix=None,ma_fix=None):
    C=np.array(d[t]['C']); r=C[:,0]; C=C[:,1]; k=kappa(float(t),2001); Ci0=C[-1]
    tail=np.abs(C-Ci0); rmax=int(r[np.where(tail>1e-11*Ci0)[0].max()]); m=(r>=rmin)&(r<=rmax); r=r[m]; C=C[m]
    best=None
    for sA in [1,-1]:
      for sB in [1,-1]:
        def f(p):
            Ci,lA,a,ms,lB,b,ma=p
            if ms_fix: ms=ms_fix*k
            if ma_fix: ma=ma_fix*k
            return (Ci+sA*np.exp(lA-a*np.log(r)-ms*r)+sB*(-1)**r*np.exp(lB-b*np.log(r)-ma*r)-C)/(np.abs(C-Ci0)+1e-12)
        s=least_squares(f,[Ci0,-2,1,2*k,-2,1,2*k],bounds=([Ci0-1e-3,-40,-3,.2*k,-40,-3,.2*k],[Ci0+1e-3,5,8,8*k,5,8,8*k]),xtol=1e-14,ftol=1e-14,max_nfev=20000)
        if best is None or s.cost<best[0]: best=(s.cost,s.x)
    return best,rmax
for t in sys.argv[1:]:
    k=kappa(float(t),2001)
    for lab,ms,ma in [('free',None,None),('ms=2k',2,None),('ms=4k',4,None),('both=2k',2,2)]:
        (c,p),rmax=run(t,8,ms,ma)
        print('t=%s %-8s cost=%.2e  ms/k=%.3f a=%.2f  ma/k=%.4f b=%.2f  (r 8..%d)'%(t,lab,c,(p[3]/k if ms is None else ms),p[2],(p[6]/k if ma is None else ma),p[5],rmax))
