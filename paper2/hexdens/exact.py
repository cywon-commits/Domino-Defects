import sys, json, numpy as np
sys.path.insert(0,'/tmp/claude-0/-home-claude/76b76c06-8fe9-5f68-b413-27ba31b71e09/scratchpad/repo/domino-defects-code/paper2')
from honey_defect import Flake
out={}
for R in [int(x) for x in sys.argv[1:]]:
    F=Flake(R); lz=np.array([F.logZ_minus(v) for v in range(len(F.A))])
    out[R]={'A':F.A.tolist(),'lz':lz.tolist()}; print(R,len(F.A),flush=True)
json.dump(out,open('exact_%s.json'%'_'.join(sys.argv[1:]),'w'))
