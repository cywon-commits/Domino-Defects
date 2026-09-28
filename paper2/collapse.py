import json, numpy as np
d=json.load(open('../data/paper2_ratio_fine.json'))
ts=[t for t in d if float(t)>0]
def cost(nu, which=1):
    curves=[]
    for t in ts:
        a=np.array(d[t]); N=a[:,0]; y=np.log(a[:,which]*N**0.25)
        curves.append((np.log(N*float(t)**nu), y))
    # pairwise overlap: interpolate each curve onto others' x in overlap region
    s=0; k=0
    for i,(xi,yi) in enumerate(curves):
        for j,(xj,yj) in enumerate(curves):
            if i==j: continue
            m=(xj>=xi.min())&(xj<=xi.max())
            if m.sum()==0: continue
            s+=np.sum((np.interp(xj[m],xi,yi)-yj[m])**2); k+=m.sum()
    return s/max(k,1), k
for which,name in [(1,'corner A (0,0)'),(2,'corner B (0,N-1)')]:
    nus=np.linspace(0.6,1.6,101); c=[cost(n,which)[0] for n in nus]
    b=nus[int(np.argmin(c))]
    print(name,'best nu=%.2f'%b,' cost(nu=1)=%.2e  cost(best)=%.2e'%(cost(1.0,which)[0],min(c)))
# asymptotic check: saturated R at small t vs t^(1/4) using largest N where saturated
d2=json.load(open('../data/paper2_ratio.json'))
for t in ['0.2','0.5','1.0']: print('t',t,'R_inf(A)=%.4f  R_inf/t^0.25=%.4f'%(d2[t][-1][1],d2[t][-1][1]/float(t)**0.25))
