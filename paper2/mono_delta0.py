import json,sys,numpy as np; sys.path.insert(0,'.'); from symbol_xi import kappa
d=json.load(open('../data/paper2_pair_big.json'))
print(' t     2kappa    delta0 (odd r)        alpha     alt. candidate would need delta0')
res={}
for t in ['0.1','0.15','0.2','0.25','0.3','0.35']:
    C=np.array(d[t]['C']); r=C[:,0].astype(int); C=C[:,1]; k=kappa(float(t),2001)
    q=np.exp(-4*k); Ci=np.mean([C[-1]-(C[-3]-C[-1])*q/(1-q), C[-2]-(C[-4]-C[-2])*q/(1-q)])
    P=C-Ci; m=(r%2==1)&(np.abs(P)>1e-10*Ci)
    rr=r[m]; Y=np.log(np.abs(P[m]))+2*k*rr
    dl=-np.diff(Y)/np.diff(rr); rm=(rr[1:]+rr[:-1])/2
    w=(rm>=rm.max()*0.3)&(rm<=rm.max()*0.9)
    X=np.c_[1/rm[w],np.ones(w.sum())]; c,resid,*_=np.linalg.lstsq(X,dl[w],rcond=None)
    # bootstrap-ish error: vary window
    alts=[]
    for lo,hi in [(0.2,0.8),(0.3,0.9),(0.4,0.95),(0.25,0.7)]:
        ww=(rm>=rm.max()*lo)&(rm<=rm.max()*hi); cc,*_=np.linalg.lstsq(np.c_[1/rm[ww],np.ones(ww.sum())],dl[ww],rcond=None); alts.append(cc[1])
    cand=2*np.arctanh(float(t)/np.sqrt(1-float(t)**2))-2*k
    res[t]=dict(delta0=float(c[1]),spread=float(np.ptp(alts)),alpha=float(c[0]),cand=float(cand))
    print(' %-5s %.4f   %+.5f ± %.5f     %.3f     %+.5f'%(t,2*k,c[1],np.ptp(alts)/2,c[0],cand))
json.dump(res,open('../data/paper2_mono_delta0.json','w'),indent=1)
