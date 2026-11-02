import sys, json, numpy as np; sys.path.insert(0,'.'); from pftri import logZ, logZ_multi
def Cinf(t): return 0.5*np.sqrt(t/(2*t*(2+t*t)+(1+2*t*t)*np.sqrt(2+t*t)))
def eta2(t):
    mu=np.sqrt(1-4*t*t+0j); x=2-mu-2*np.sqrt(1-t*t-mu); return np.real(-0.5*np.log(x))
L=64; out={}
for t in [0.2,0.3,0.4,0.5]:
    lz=logZ(L,L,t); c0=L//2; res={}
    for name,step in [('row',(0,1)),('diag',(1,1)),('anti',(1,-1))]:
        ks=range(1,13) if name!='row' else range(1,25); d=[]
        for k in ks:
            dr,dc=step[0]*k,step[1]*k; u=(c0-dr//2,c0-dc//2); v=(u[0]+dr,u[1]+dc)
            d.append(np.exp(logZ_multi(L,L,t,[u,v])-lz)-Cinf(t))
        d=np.array(d); man=np.array([abs(step[0])*k+abs(step[1])*k for k in ks])
        res[name]=dict(man=man.tolist(),d=d.tolist())
        m=man>=8; a=np.abs(d)
        sl=np.polyfit(man[m],np.log(a[m]),1)[0] if m.sum()>2 else np.nan
        res[name]['rate_per_manhattan']=-sl
    out[t]=dict(res=res,two_ln_eta2=2*eta2(t))
    print(t,'2ln eta2=%.3f'%(2*eta2(t)),{k:round(v['rate_per_manhattan'],3) for k,v in res.items()},flush=True)
json.dump({str(k):v for k,v in out.items()},open('../data/paper2_antidiag.json','w'))
