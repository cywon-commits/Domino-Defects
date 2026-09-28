import sys,json,time,numpy as np; sys.path.insert(0,'.'); from pairdirect import PairD
out={}; T0=time.time()
for t,L in [(0.15,256),(0.2,256),(0.25,256),(0.3,256),(0.35,256),(0.4,256),(0.45,256),(0.5,256),(0.6,256),(0.7,256),(1.0,256),(0.1,384)]:
    P=PairD(L,t); row=[]
    for r in range(1,81):
        c1=L//2-(r+1)//2; row.append([r,P.C(c1,c1+r)])
    out[str(t)]={'L':L,'C':row}; json.dump(out,open('../data/paper2_pair_big.json','w'))
    print('t=%.2f L=%d C(80)=%.14f  %.0fs'%(t,L,row[-1][1],time.time()-T0),flush=True)
