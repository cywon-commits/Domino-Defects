import numpy as np, sys, scipy.sparse as sp, scipy.sparse.linalg as spl
sys.path.insert(0,'.'); from honeyflake import flake; from honey_meanh import meanh
from scipy.interpolate import RegularGridInterpolator
# u: harmonic in the regular hexagon (apothem 1, vertices at 30+60k), u = +1 on edges with normal 0,120,240, -1 on 60,180,300
n=401; x=np.linspace(-1.2,1.2,n); hgrid=x[1]-x[0]; XX,YY=np.meshgrid(x,x,indexing='ij')
nrm=[np.array([np.cos(np.radians(60*k)),np.sin(np.radians(60*k))]) for k in range(6)]
proj=np.stack([XX*v[0]+YY*v[1] for v in nrm]); inside=np.all(proj<1,axis=0)
edge=np.argmax(proj,axis=0); bval=np.where(edge%2==0,1.0,-1.0)
idx=-np.ones((n,n),int); idx[inside]=np.arange(inside.sum()); N=inside.sum()
I=[];J=[];V=[];b=np.zeros(N)
for i,j in zip(*np.where(inside)):
    r=idx[i,j]; I.append(r);J.append(r);V.append(-4.0)
    for di,dj in [(1,0),(-1,0),(0,1),(0,-1)]:
        ii,jj=i+di,j+dj
        if inside[ii,jj]: I.append(r);J.append(idx[ii,jj]);V.append(1.0)
        else: b[r]-=bval[ii,jj]
U=np.zeros((n,n)); U[inside]=spl.spsolve(sp.csr_matrix((V,(I,J)),shape=(N,N)),b); U[~inside]=bval[~inside]
u=RegularGridInterpolator((x,x),U)
fams={'J=0':[12,17,22,25,27,30],'J=+-1':[11,16,21,24,26,29],'J=+-2':[15,18,20,23,28]}
for fam,Rs in fams.items():
    for R in Rs:
        A,B=flake(R,'armchair')
        if len(A)!=len(B): continue
        P,H,err=meanh(A,B); c0=np.array([0,-1.0]); Z=(P-c0)/R
        rr=np.array([max(z@v for v in nrm) for z in Z])        # "hexagonal radius" (1 at the boundary)
        out=[]
        for rmax in [0.5,0.7,0.85]:
            m=rr<rmax; X=np.c_[np.ones(m.sum()),u(Z[m])]; c,res,*_=np.linalg.lstsq(X,H[m],rcond=None)
            out.append(c[1])
        a=out[1]
        # J at the 30-degree corner = level(edge 60) - level(edge 0) = -a - a
        print('%-6s R=%2d  a(rmax=.5,.7,.85)= %s   J_eff(30deg corner) = %+.3f'%(fam,R,' '.join('%+.3f'%v for v in out),-2*a),flush=True)
