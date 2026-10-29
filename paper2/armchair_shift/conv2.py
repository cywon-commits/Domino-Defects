import sys, json, numpy as np; sys.path.insert(0,'.'); import delta as dl
from armstrip import build, probs
out=[]
for R in [float(x) for x in sys.argv[1:]]:
    A,B,bonds=build(R); P1,_=probs(A,B,bonds,1024); P2,_=probs(A,B,bonds,2048); c=np.abs(P1-P2).max()
    row=[R,len(A),c]+[dl.delta(R,dcut=f*R,pre=(A,B,bonds,P2))['diff'] for f in (0.15,0.3,0.45)]
    out.append(row); print('R=%5g nA=%4d k-conv=%.0e  diff = %+.6f %+.6f %+.6f'%tuple(row),flush=True)
json.dump(out,open('conv2.json','w'))
