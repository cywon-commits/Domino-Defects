import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import numpy as np, math
from var2 import Mesh
from kreg import aztec
n=32; M=Mesh(aztec(n)); pure0,h0,exc,r=M.solve(schedule=((50,1e-3),(500,1e-4)),iters=20000)
g=M.grads(h0); P=M.pts; C=(P[M.tri[:,0]]+P[M.tri[:,1]]+P[M.tri[:,2]])/3
xh=(C[:,0]-n)/n; yh=(C[:,1]-n)/n
def cep(x,y):
    e=np.sqrt(np.maximum(1-2*x*x-2*y*y,1e-12)); f=lambda a: 0.5+np.arctan(a/e)/np.pi
    return f(2*x-1),f(-2*x-1),f(2*y-1),f(-2*y-1)
sel=(xh**2+yh**2<0.3)
E,W,N,S=cep(xh[sel],yh[sel])
for name,(a,b) in {'s vs E-W':(g[sel,0],E-W),'t vs N-S':(g[sel,1],N-S),'s vs N-S':(g[sel,0],N-S),'t vs E-W':(g[sel,1],E-W)}.items():
    print(name, 'corr=%.3f'%np.corrcoef(a,b)[0,1], 'slope=%.3f'%np.polyfit(b,a,1)[0])
