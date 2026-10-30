import json, numpy as np, sys; sys.path.insert(0,'.'); from sc import density, inv
d=json.load(open('exact_28_27_36_38.json'))
rules={'I':{'edge-class':[-1/6]*6,'naive':[0,-1/3]*3},
       'II':{'edge-class':[1/6,-1/2]*3,'naive':[1/3,-2/3]*3}}
typ={28:'I',36:'I',27:'II',38:'II'}
res={}
for R in [28,36,27,38]:
    A=np.array(d[str(R)]['A']); lz=np.array(d[str(R)]['lz']); Z=A[:,0]+1j*A[:,1]
    hexr=np.max(np.stack([A@np.array([np.cos(np.radians(60*k)),np.sin(np.radians(60*k))]) for k in range(6)]),axis=0)
    row=[]
    for c in [0.0,0.5,1.0]:
        a=R+c; pre=inv(Z,a)
        for rule,q in rules[typ[R]].items():
            lf,_=density(Z,a,q,pre=pre)
            for depth in [3,6]:
                m=hexr<R-depth
                diff=lz[m]-lf[m]; diff-=diff.mean()
                row.append((c,rule,depth,np.sqrt(np.mean(diff**2)),np.max(np.abs(diff))))
    res[R]=row
    print(flush=True); print('R=%d type %s  (rms / max deviation of ln rho, interior cut at depth 3 and 6)'%(R,typ[R]))
    for c,rule,dep,rms,mx in row: print('   a=R+%.1f  %-10s depth %d: rms %.3f  max %.3f'%(c,rule,dep,rms,mx))
json.dump({str(k):v for k,v in res.items()},open('compare.json','w'))
