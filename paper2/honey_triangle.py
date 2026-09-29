"""60-degree armchair corners (the Ciucu-Fischer boundary type): C3-symmetric armchair triangles
centred on an A site. Forced defects; corner exponent vs the symmetry-forced prediction."""
import numpy as np, sys, json; sys.path.insert(0,'.')
import honey_defect as hd
from honey_corner import boundary
s3=np.sqrt(3); D=hd.D
def triA(R):
    a1=np.array([s3,0]); a2=np.array([s3/2,1.5]); M=int(2.2*R)+4
    A=np.array([m*a1+n*a2 for m in range(-M,M+1) for n in range(-M,M+1)]); B=A+D[0]
    nrm=[np.array([np.cos(np.radians(a)),np.sin(np.radians(a))]) for a in (0,120,240)]
    ins=lambda P: np.all([P@v<=R for v in nrm],axis=0); A=A[ins(A)]; B=B[ins(B)]
    key=lambda p:(round(p[0],4),round(p[1],4))
    while True:
        Bs={key(b) for b in B}; As={key(a) for a in A}
        dA=np.array([sum(key(a+d) in Bs for d in D) for a in A]); dB=np.array([sum(key(b-d) in As for d in D) for b in B])
        if (dA>=2).all() and (dB>=2).all(): break
        A=A[dA>=2]; B=B[dB>=2]
    return A,B
class Tri(hd.Flake):
    def __init__(s,R): s.A,s.B=triA(R); s.E=hd.bonds(s.A,s.B); s.R=R
if __name__=='__main__':
    ok=[]
    for R in np.arange(6,60,1.0):
        A,B=triA(R); d=len(A)-len(B)
        if d==1:
            X,H,tot=boundary(A,B); ang=np.degrees(np.arctan2(X[:,1],X[:,0]))%360
            # edges have normals 0,120,240; corners at 60,180,300.  average over the middle of each edge
            lev=[]
            for nang in (0,120,240):
                dd=np.abs(((ang-nang+180)%360)-180); lev.append(H[dd<25].mean())
            J=[lev[1]-lev[0],lev[2]-lev[1],lev[0]-lev[2]]
            ok.append(R); print('R=%4.1f nA=%d nA-nB=%d winding=%d  naive J at corners 60,180,300: %s'%(R,len(A),d,tot,' '.join('%+.2f'%j for j in J)))
    json.dump(ok,open('/tmp/triR.json','w'))
