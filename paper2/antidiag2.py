import sys, json, numpy as np; sys.path.insert(0,'.'); from pftri import logZ, logZ_multi
def Cinf(t): return 0.5*np.sqrt(t/(2*t*(2+t*t)+(1+2*t*t)*np.sqrt(2+t*t)))
out={}
for t,L in [(0.5,96),(0.3,96)]:
    lz=logZ(L,L,t); c0=L//2; d=[]
    for k in range(1,11):
        u=(c0-k//2,c0+k//2); v=(u[0]+k,u[1]-k)
        d.append(np.exp(logZ_multi(L,L,t,[u,v])-lz)-Cinf(t))
    out[str(t)]=d; print(t,L,['%.3e'%x for x in d],flush=True)
json.dump(out,open('../data/paper2_antidiag96.json','w'))
