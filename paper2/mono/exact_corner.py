import sys,json,numpy as np; sys.path.insert(0,'../trapping'); from trap import Board
out={}
for L in [17,25,33,49,65,97,129]:
    B=Board(L); f=lambda r,c: np.exp(B.logZ([(r,c)]))
    corner=[(r,c) for r in range(3) for c in range(3) if r+c<=2 and (r+c)%2==0]
    h=L//2; cen=[(h+a,h+b) for a in range(-2,3) for b in range(-2,3) if abs(a)+abs(b)<=2 and (a+b)%2==0]
    lc=np.array([B.logZ([p]) for p in corner]); le=np.array([B.logZ([p]) for p in cen]); m=max(lc.max(),le.max())
    rat=np.mean(np.exp(lc-m))/np.mean(np.exp(le-m)); rat1=np.exp(B.logZ([(0,0)])-B.logZ([(h,h)]))
    out[L]=dict(avg=rat,single=rat1); print(L,rat,rat1,flush=True)
json.dump(out,open('../../data/paper2_exact_corner_bv0.json','w'))
Ls=np.array(sorted(out)); a=np.array([out[L]['avg'] for L in Ls]); s=np.array([out[L]['single'] for L in Ls])
for name,y in [('avg',a),('single',s)]:
    m=(Ls>=25)&(Ls<=65); print(name,'gamma(25-65)=%.4f'%-np.polyfit(np.log(Ls[m]),np.log(y[m]),1)[0], 'local',[round(-np.log(y[i+1]/y[i])/np.log(Ls[i+1]/Ls[i]),4) for i in range(len(Ls)-1)])
