import numpy as np, sys, json; sys.path.insert(0,'.')
import honey_defect as hd; from honey_corner import boundary
from honey_para import paraA, n0, n60, corners, SHIFT, PIOT
def exactJ(A,B,c0):
    X,H,tot=boundary(A,B); P=X-c0; n=len(H)
    cp=corners[(-1,1)][0]; i0=int(np.argmin(np.linalg.norm(P-cp*np.abs(P@n0).max(),axis=1)))  # start at corner 120
    k=np.arange(n); H2=H[(i0+k)%n]+np.where(i0+k>=n,tot,0); P2=P[(i0+k)%n]
    lev={}
    for lab,nv,sgn in [('0',n0,1),('60',n60,1),('180',n0,-1),('240',n60,-1)]:
        s=sgn*(P2@nv); other=P2@(n60 if lab in ('0','180') else n0)
        m=(s>0.9*s.max())&(np.abs(other)<0.5*np.abs(other).max()); lev[lab]=np.median(H2[m])
    J={30:lev['60']-lev['0'],120:lev['180']+tot-lev['60'],210:lev['240']-lev['180'],300:lev['0']-lev['240']}
    return J,tot
class Para(hd.Flake):
    def __init__(s,R): s.A,s.B=paraA(R,np.zeros(2)); s.E=hd.bonds(s.A,s.B); s.R=R
fams={'F1':[9,12,14,17,22,25,27,30,35,38],'F2':[10,15,18,20,23,28,31,33,36]}
res={}
for fam,Rs in fams.items():
    rows=[]
    for R in Rs:
        F=Para(R); A=F.A; J,tot=exactJ(F.A,F.B,np.zeros(2))
        lz0=F.logZ_minus(int(np.argmin(np.linalg.norm(A,axis=1)))); r=[]
        for c,(p,ang) in [(c,corners[k]) for c,k in [(30,(1,1)),(120,(-1,1)),(210,(-1,-1)),(300,(1,-1))]]:
            u=p/np.linalg.norm(p); s=A@u; perp=np.linalg.norm(A-np.outer(s,u),axis=1); m=np.where(perp<0.9)[0]
            r.append(np.exp(F.logZ_minus(m[np.argmax(s[m])])-lz0))
        rows.append([R]+r); Jr={c:round(v) for c,v in J.items()}
    a=np.array(rows)
    e=+1; pred={}
    for c in (30,120,210,300):
        q=-(Jr[c]+SHIFT[c])/3; qn=-Jr[c]/3
        pred[c]=(PIOT[c]*(0.5+e*q)-0.25, PIOT[c]*(0.5+e*qn)-0.25)
    print('%s: exact naive J (sum %d) = %s'%(fam,sum(Jr.values()),Jr))
    for j,c in enumerate((30,120,210,300)):
        y=a[:,j+1]; g=-np.polyfit(np.log(a[:,0]),np.log(y),1)[0]; g2=-np.polyfit(np.log(a[-5:,0]),np.log(y[-5:]),1)[0]
        print('   corner %3d (%s): predicted gamma %+.3f  [naive rule %+.3f]   measured %+.3f (all R), %+.3f (last 5 R)'%(c,'120deg' if PIOT[c]==1.5 else ' 60deg',pred[c][0],pred[c][1],g,g2),flush=True)
    res[fam]={'J':Jr,'data':a.tolist()}
json.dump(res,open('../data/paper2_honey_para.json','w'))
