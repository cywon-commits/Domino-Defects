import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np, math, json, time
from thmB import Quadrant, ratio
Q=Quadrant(6000)
res={'diag':[], 'ang':[]}
t=time.time()
for r in [2,4,8,16,32,64,128,256]:
    v=(r,r+1); R=ratio(Q.G,v); rho=math.hypot(r+0.5,r+1.5)
    res['diag'].append([r,rho,R]); print("diag",r,rho,R,flush=True)
prev=None
for r,rho,R in res['diag']:
    if prev: print("  local exponent %.4f"%(math.log(R/prev[1])/math.log(rho/prev[0])))
    prev=(rho,R)
# angular profile at radius ~120
Rr=120
for ang in np.linspace(5,85,17):
    a=math.radians(ang); x=int(round(Rr*math.cos(a)-0.5)); y=int(round(Rr*math.sin(a)-0.5))
    if (x+y)%2==0: x+=1
    R=ratio(Q.G,(y,x)); rho=math.hypot(x+0.5,y+0.5); phi=math.atan2(y+0.5,x+0.5)
    res['ang'].append([phi,rho,R]); print("ang",round(math.degrees(phi),1),rho,R,flush=True)
json.dump(res,open("thmB_res.json","w")); print("time",time.time()-t)
