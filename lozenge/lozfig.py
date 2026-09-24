import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import sys, json, math
sys.argv=['x']; exec(open('rate.py').read().split("cases=")[0])
import numpy as np
out={}
for al,be,ga in [(1,1,1),(1,1.5,0.75)]:
    cc=dlogZc(al,be,ga); rows=[]
    for kap in [0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5]:
        lam=kap
        if lam>=be or kap>=al: continue
        p1,pm,tb=predict(al,be,ga,kap,lam,cc)
        vals=[]
        for n in (60,120,240):
            a,b,c=int(al*n),int(be*n),int(ga*n); k=int(kap*n)+1; l=int(lam*n)+1
            vals.append(float(mp.log(lgvW(a,b,c,k,l)/lgvZ(a,b,c))/n))
        ns=np.array([60,120,240.]); A=np.stack([np.ones(3),np.log(ns)/ns,1/ns],1); ex=float(np.linalg.solve(A,np.array(vals))[0])
        rows.append(dict(kap=kap,tau1=p1,rule=pm,tau=float(tb),vals=vals,extrap=ex)); print(al,be,ga,rows[-1],flush=True)
    out[f'{al},{be},{ga}']=dict(cc=cc,rows=rows)
json.dump(out,open('lozfig.json','w'))
