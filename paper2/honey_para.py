"""Prediction test of the edge-class shift: armchair parallelograms (two 60 and two 120 degree corners).
Edges with normals 0 and 240 (class +) and 60 and 180 (class -)."""
import numpy as np, sys, json; sys.path.insert(0,'.')
import honey_defect as hd; from honey_corner import boundary
s3=np.sqrt(3); D=hd.D
n0=np.array([1.0,0]); n60=np.array([0.5,s3/2])
def paraA(R,c0):
    a1=np.array([s3,0]); a2=np.array([s3/2,1.5]); M=int(2.6*R)+4
    A=np.array([m*a1+n*a2 for m in range(-M,M+1) for n in range(-M,M+1)]); B=A+D[0]
    ins=lambda P:(np.abs((P-c0)@n0)<=R)&(np.abs((P-c0)@n60)<=R); A=A[ins(A)]; B=B[ins(B)]
    key=lambda p:(round(p[0],4),round(p[1],4))
    while True:
        Bs={key(b) for b in B}; As={key(a) for a in A}
        dA=np.array([sum(key(a+d) in Bs for d in D) for a in A]); dB=np.array([sum(key(b-d) in As for d in D) for b in B])
        if (dA>=2).all() and (dB>=2).all(): break
        A=A[dA>=2]; B=B[dB>=2]
    return A,B
# corners: P.n0 = s1 R, P.n60 = s2 R
Mx=np.array([n0,n60]); corners={}
for s1,s2 in [(1,1),(-1,1),(-1,-1),(1,-1)]:
    p=np.linalg.solve(Mx,np.array([s1,s2]))
    ang=(np.degrees(np.arctan2(p[1],p[0]))%360)
    interior=60 if s1==s2 else 120      # (+,+): between normals 0 and 60 -> 120? check below
    corners[(s1,s2)]=(p,ang)
if __name__=='__main__':
    for (s1,s2),(p,ang) in corners.items(): print('corner',(s1,s2),'dir %.1f'%ang, 'unit pos',np.round(p,3))

def naiveJ(A,B,c0):
    X,H,tot=boundary(A,B); P=X-c0
    lev={}
    for lab,n,sgn in [('0',n0,1),('60',n60,1),('180',n0,-1),('240',n60,-1)]:
        s=sgn*(P@n); other=P@(n60 if lab in ('0','180') else n0)
        m=(s>0.9*s.max())&(np.abs(other)<0.5*np.abs(other).max()); lev[lab]=np.median(H[m])
    # ccw order of edges: 0 -> [corner 30] -> 60 -> [corner 120] -> 180 -> [corner 210] -> 240 -> [corner 300] -> 0
    J={30:lev['60']-lev['0'],120:lev['180']-lev['60'],210:lev['240']-lev['180'],300:lev['0']-lev['240']}
    return J,tot
SHIFT={30:+0.5,120:0.0,210:-0.5,300:0.0}; PIOT={30:1.5,120:3.0,210:1.5,300:3.0}
def predict(J,e):
    out={}
    for c,j in J.items():
        jj=j
        # resolve the mod-3 ambiguity from the winding cut: pick the representative in [-1.5,1.5]+... closest to |j|<=2
        while jj>2.25: jj-=3
        while jj<-2.25: jj+=3
        q=-(jj+SHIFT[c])/3; out[c]=(jj,q,PIOT[c]*(0.5+e*q)-0.25)
    return out
if __name__=='__main__' and len(sys.argv)>1:
    fam={}
    for R in np.arange(8,46,1.0):
        for cname,c0 in [('A',np.zeros(2)),('face',np.array([0,-1.0]))]:
            A,B=paraA(R,c0); d=len(A)-len(B)
            if abs(d)!=1: continue
            J,tot=naiveJ(A,B,c0); e=+1 if d==1 else -1
            pr=predict(J,e); sig=(cname,e)+tuple(round(pr[c][0]) for c in (30,120,210,300))
            fam.setdefault(sig,[]).append(R)
    for sig,Rs in fam.items():
        print('centre=%s e=%+d naive J(30,120,210,300)=%s  R=%s'%(sig[0],sig[1],sig[2:],Rs[:12]))
    json.dump({str(k):v for k,v in fam.items()},open('/tmp/parafam.json','w'))
