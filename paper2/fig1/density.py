import numpy as np, json, sys; sys.path.insert(0,'.'); from trap import Board
L=24; B=Board(L); W=[(r,c) for r in range(L) for c in range(L) if (r+c)%2]
lz=np.array([B.logZ([(0,0),v]) for v in W]); p=np.exp(lz-lz.max()); p/=p.sum()
json.dump({'L':L,'W':W,'p':p.tolist()},open('density.json','w')); print('max/min',p.max()/p.min())
