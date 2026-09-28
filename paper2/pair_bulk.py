"""Task 1 step 2: bulk pair correlation C(r)=Z(board-{u,v})/Z(board) on L x L, pair on middle row,
placed symmetrically about the center column."""
import sys, json, time, numpy as np; sys.path.insert(0,'.'); from pftri import logZ, logZ_multi
L=64; r0=L//2; ts=[0,0.1,0.2,0.3,0.5,0.7,1.0]; rs=list(range(1,33))
out={}; T0=time.time()
for t in ts:
    lz=logZ(L,L,t); row=[]
    for r in rs:
        c1=L//2-(r+1)//2; c2=c1+r
        row.append([r, float(np.exp(logZ_multi(L,L,t,[(r0,c1),(r0,c2)])-lz))])
    out[str(t)]=row
    print(t, ' '.join('%d:%.4g'%(a,b) for a,b in row[::4]), '%.0fs'%(time.time()-T0), flush=True)
json.dump(out,open('../data/paper2_pair_bulk.json','w'),indent=1)
