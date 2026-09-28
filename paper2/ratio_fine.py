"""R(N,t) for N = 1 mod 4, fine t grid near 0. Corners: A=(0,0) (diagonal incident), B=(0,N-1)."""
import sys, json, time, numpy as np; sys.path.insert(0,'.'); from pftri import logZ
ts=[0,0.005,0.01,0.02,0.03,0.05,0.07,0.1,0.15,0.2]
Ns=[13,21,29,37,45,53,61,69,81]
out={}; t0=time.time()
for t in ts:
    rows=[]
    for N in Ns:
        m=(N-1)//2; lm=logZ(N,N,t,(m,m))
        rows.append([N,float(np.exp(logZ(N,N,t,(0,0))-lm)),float(np.exp(logZ(N,N,t,(0,N-1))-lm))])
    out[str(t)]=rows
    print(t,' '.join('%d:%.4f'%(r[0],r[1]*r[0]**0.25) for r in rows),'%.0fs'%(time.time()-t0),flush=True)
json.dump(out,open('../data/paper2_ratio_fine.json','w'),indent=1)
