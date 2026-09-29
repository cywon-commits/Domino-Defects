import numpy as np, sys, json; sys.path.insert(0,'.'); from honey_triangle import Tri
rows=[]
for R in [8,12,16,20,24,28,32,36,40,46,52,58]:
    F=Tri(R); A=F.A; c=int(np.argmin(np.linalg.norm(A,axis=1))); lz0=F.logZ_minus(c)
    vals=[]
    for ang in (60,180,300):
        u=np.array([np.cos(np.radians(ang)),np.sin(np.radians(ang))]); s=A@u
        perp=np.linalg.norm(A-np.outer(s,u),axis=1); m=np.where(perp<0.9)[0]; m=m[np.argsort(-s[m])]
        vals.append([np.exp(F.logZ_minus(i)-lz0) for i in m[:3]])
    rows.append([R,len(A)]+[v[0] for v in vals]+vals[0][1:])
    print('R=%d nA=%d  corner/centre (60,180,300) = %s   next sites on the 60 line: %s'%(R,len(A),' '.join('%.4g'%v[0] for v in vals),' '.join('%.4g'%x for x in vals[0][1:])),flush=True)
a=np.array(rows); json.dump(a.tolist(),open('../data/paper2_honey_triangle.json','w'))
for j,lab in [(2,'corner site'),(5,'2nd site'),(6,'3rd site')]:
    p=np.polyfit(np.log(a[:,0]),np.log(a[:,j]),1); loc=np.diff(np.log(a[:,j]))/np.diff(np.log(a[:,0]))
    print('%-11s gamma = %.3f   local: %s'%(lab,-p[0],' '.join('%.3f'%-x for x in loc)))
