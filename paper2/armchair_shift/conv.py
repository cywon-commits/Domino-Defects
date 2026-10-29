import sys, json, numpy as np; sys.path.insert(0,'.'); import delta as dl
from armstrip import build, probs
out=[]
for R in [float(x) for x in sys.argv[1:]]:
    A,B,bonds=build(R)
    if len(A)!=len(B): print('R=%g unbalanced'%R); continue
    P1,_=probs(A,B,bonds,192); P2,_=probs(A,B,bonds,384); c=np.abs(P1-P2).max()
    row=[R,len(A),c]
    for frac in [0.15,0.3,0.45]:
        r=dl.delta(R,dcut=frac*R,pre=(A,B,bonds,P2)); row.append(r['diff'])
    out.append(row); print('R=%5g nA=%4d k-conv=%.0e  diff(cut 0.15R,0.3R,0.45R) = %+.6f %+.6f %+.6f'%tuple(row),flush=True)
    json.dump(out,open('conv_%s.json'%sys.argv[1],'w'))
