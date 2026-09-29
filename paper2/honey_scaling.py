import numpy as np, sys, json; sys.path.insert(0,'.'); from honey_defect import Flake
def corner_site(F,R,cang):
    u=np.array([np.cos(np.radians(cang)),np.sin(np.radians(cang))]); s=F.A@u
    perp=np.linalg.norm(F.A-np.outer(s,u),axis=1); m=np.where(perp<0.9)[0]; return m[np.argmax(s[m])]
res={}
for typ,Rs in [('I',[10,15,18,20,23,28,31,33,36]),('II',[9,12,14,17,22,25,27,30,35,38])]:
    rows=[]
    for R in Rs:
        F=Flake(R); lz0=F.logZ_minus(int(np.argmin(np.linalg.norm(F.A,axis=1))))
        r30=np.exp(F.logZ_minus(corner_site(F,R,30))-lz0); r90=np.exp(F.logZ_minus(corner_site(F,R,90))-lz0)
        rows.append((R,r30,r90))
    a=np.array(rows); res[typ]=a.tolist()
    for j,c in [(1,30),(2,90)]:
        p=np.polyfit(np.log(a[:,0]),np.log(a[:,j]),1); loc=np.diff(np.log(a[:,j]))/np.diff(np.log(a[:,0]))
        print('type %-2s corner %d:  ratio(R) %s   fitted gamma = %.3f   local gammas(last 3) %s'%(typ,c,' '.join('%.3g'%x for x in a[:,j]),-p[0],' '.join('%.3f'%-x for x in loc[-3:])),flush=True)
json.dump(res,open('../data/paper2_honey_scaling.json','w'))
