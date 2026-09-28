import sys, json, time, numpy as np; sys.path.insert(0,'.'); from pairfast import Pair
L=96; out={}; T0=time.time()
for t in [0.1,0.2,0.3,0.5,0.7,1.0]:
    P=Pair(L,t); row=[]
    for r in range(1,49):
        c1=L//2-(r+1)//2; row.append([r,float(P.C(c1,c1+r))])
    out[str(t)]=row; del P
    json.dump(out,open('../data/paper2_pair96.json','w'))
    print(t,'C(48)=%.12f'%row[-1][1],'%.0fs'%(time.time()-T0),flush=True)
