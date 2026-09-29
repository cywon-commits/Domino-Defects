import numpy as np, sys, json, time; sys.path.insert(0,'.'); from honey_defect import Flake
from honey_corner import boundary
def jumps(A,B):
    X,H,tot=boundary(A,B); ang=np.degrees(np.arctan2(X[:,1],X[:,0]))%360
    # put the branch cut at the 150-degree corner: re-base heights so the cycle starts there
    i0=np.argmin(np.abs(((ang-150+180)%360)-180)); H=np.roll(H,-i0)-H[i0]; ang=np.roll(ang,-i0)
    H=np.where(np.arange(len(H))>0,H,H)  # heights now continuous away from 150 deg
    a2=np.where(ang<30,ang+360,ang); avg=[]
    for k in range(6):
        lo=30+60*k; m=(a2>lo+12)&(a2<lo+48); avg.append(H[m].mean())
    return [avg[k]-avg[k-1] for k in range(6)],tot
out={}; T0=time.time()
Rs=[int(x) for x in sys.argv[1:]]
for R in Rs:
    F=Flake(R); A=F.A
    J,tot=jumps(F.A,F.B)
    lz0=F.logZ_minus(int(np.argmin(np.linalg.norm(A,axis=1))))      # centre (A at origin)
    rec={'J':J,'corner':[],'edge':[]}
    for k in range(6):
        ca=np.radians(30+60*k); cp=R/np.cos(np.pi/6)*np.array([np.cos(ca),np.sin(ca)])
        ea=np.radians(60*k);   ep=R*np.array([np.cos(ea),np.sin(ea)])
        for lab,pt,rad in [('corner',cp,2.6),('edge',ep,1.8)]:
            idx=np.where(np.linalg.norm(A-pt,axis=1)<rad)[0]
            s=sum(np.exp(F.logZ_minus(v)-lz0) for v in idx)
            rec[lab].append(float(s))
    out[R]=rec
    print('R=%d  J=%s  corner/centre=%s  edge/centre=%s  (%.0fs)'%(R,' '.join('%+.1f'%j for j in J),' '.join('%.4g'%x for x in rec['corner']),' '.join('%.4g'%x for x in rec['edge']),time.time()-T0),flush=True)
    json.dump(out,open('../data/paper2_honey_task2_%s.json'%('_'.join(sys.argv[1:3])),'w'))
