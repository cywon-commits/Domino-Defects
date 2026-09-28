import numpy as np, sys, json; sys.path.insert(0,'.'); from honeyflake import flake, probs
out={}
for kind in ['armchair','zigzag']:
    for R in ([12,20,28] if kind=="armchair" else [12,20]):
        A,B=flake(R,kind); E,P=probs(A,B)
        pos=np.array([(A[e[0]]+B[e[1]])/2 for e in E]); r=np.linalg.norm(pos-np.array([0,-1.0]),axis=1)/R
        dev=np.abs(P-1/3); s='%-8s R=%2d sites=%4d '%(kind,R,2*len(A))
        for lo,hi in [(0,0.5),(0.5,0.8),(0.8,0.95),(0.95,1.05)]:
            m=(r>=lo)&(r<hi); s+=' [%.2f,%.2f): mean|P-1/3|=%.3f max=%.3f'%(lo,hi,dev[m].mean(),dev[m].max())
        s+='  | edges with P<0.05 or >0.95: %d'%np.sum((P<0.05)|(P>0.95))
        print(s,flush=True); out['%s_%d'%(kind,R)]=[pos.tolist(),P.tolist()]
json.dump(out,open('../data/paper2_honeyflake.json','w'))
