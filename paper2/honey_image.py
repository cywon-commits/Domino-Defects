"""Same-sign image test: defect probability at edge-midpoint sites vs R, and inward profile."""
import numpy as np, sys, json; sys.path.insert(0,'.'); from honey_defect import Flake
def line_sites(F,ang_deg,maxd):
    n=np.array([np.cos(np.radians(ang_deg)),np.sin(np.radians(ang_deg))]); t=np.array([-n[1],n[0]])
    s=F.A@n; perp=np.abs(F.A@t); m=np.where(perp<0.9)[0]; m=m[np.argsort(-s[m])]
    return [(F.R-s[i],i) for i in m if F.R-s[i]<maxd]
res={}
for typ,Rs in [('I',[10,15,18,20,23,28,31,33,36]),('II',[9,12,14,17,22,25,27,30,35,38])]:
    rows=[]
    for R in Rs:
        F=Flake(R); lz0=F.logZ_minus(int(np.argmin(np.linalg.norm(F.A,axis=1))))
        vals=[]
        for ang in [0,60]:
            d,i=line_sites(F,ang,3)[0]; vals.append(np.exp(F.logZ_minus(i)-lz0))
        rows.append([R]+vals)
    a=np.array(rows); res[typ]=a.tolist()
    for j,ang in [(1,0),(2,60)]:
        p=np.polyfit(np.log(a[:,0]),np.log(a[:,j]),1); loc=np.diff(np.log(a[:,j]))/np.diff(np.log(a[:,0]))
        print('type %-2s edge %2d deg: edge/centre %s  exponent %.3f  local(last 4) %s'%(typ,ang,' '.join('%.3g'%x for x in a[:,j]),p[0],' '.join('%.3f'%x for x in loc[-4:])),flush=True)
# inward profile at the largest R of each type
for typ,R in [('I',36),('II',38)]:
    F=Flake(R); lz0=F.logZ_minus(int(np.argmin(np.linalg.norm(F.A,axis=1))))
    for ang in [0,60]:
        pr=[(d,np.exp(F.logZ_minus(i)-lz0)) for d,i in line_sites(F,ang,R*0.6)]
        d=np.array([x[0] for x in pr]); r=np.array([x[1] for x in pr])
        m=(d>2)&(d<R/4); sl=np.polyfit(np.log(d[m]),np.log(r[m]),1)[0] if m.sum()>3 else np.nan
        res['profile_%s_%d'%(typ,ang)]=[d.tolist(),r.tolist()]
        print('type %-2s R=%d edge %2d deg inward: slope dlnρ/dlnd over 2<d<R/4 = %.3f   first values %s'%(typ,R,ang,sl,' '.join('%.1f:%.3f'%(x,y) for x,y in pr[:6])),flush=True)
json.dump(res,open('../data/paper2_honey_image.json','w'))
