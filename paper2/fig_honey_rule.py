import json,numpy as np,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
C1,C2,GR='#1f3a5f','#dd6b20','#8a8f96'
rows=json.load(open('../data/paper2_corner_exponents_err.json'))
fig,ax=plt.subplots(1,2,figsize=(9,3.8))
seen=set()
for name,pr,nv,me,er in rows:
    if name.startswith('para F2') and name!='para F2 30': continue
    temp=name.startswith('hex I ') or name.startswith('tri') or name.startswith('para F2')
    col=C1 if temp else C2; lab=('Temperleyan (hex I, triangle, F2)' if temp else 'not Temperleyan (hex II, F1)')
    ax[0].errorbar(pr,me,yerr=er,fmt='o' if temp else 's',color=col,mfc=col if temp else 'none',ms=6,capsize=2,label=None if lab in seen else lab); seen.add(lab)
    if abs(nv-pr)>1e-9: ax[0].plot(pr,nv,'x',color=GR,ms=7,label=None if 'naive' in seen else 'naive rule at the same corner'); seen.add('naive')
xx=np.array([-0.7,1.2]); ax[0].plot(xx,xx,':',color='k',lw=0.8)
ax[0].set_xlabel(r'predicted $\gamma$, Eq. (5)'); ax[0].set_ylabel(r'measured $\gamma$'); ax[0].legend(frameon=False,fontsize=8,loc='upper left'); ax[0].set_title('(a) corner exponents',fontsize=10)
d=json.load(open('../data/paper2_armchair_strip_delta.json'))
R=np.array(sorted(int(k) for k in d)); m=np.array([-np.mean(d[str(r)]) for r in R]); s=np.array([np.ptp(d[str(r)])/2 for r in R])
ax[1].errorbar(1/R,m,yerr=s,fmt='o',color=C1,ms=5,capsize=2); ax[1].axhline(0.5,color=C2,ls='--',lw=1,label=r'$\delta_--\delta_+=1/2$')
for r,y in zip(R,m):
    if r in (12,40,320): ax[1].annotate('R = %d'%r,(1/r,y),textcoords='offset points',xytext=(-40,-3) if r==12 else (-12,-14),fontsize=8)
ax[1].set_xlabel('1/R  (R = strip half-width)'); ax[1].set_ylabel(r'$\delta_--\delta_+$'); ax[1].legend(frameon=False,fontsize=8,loc='center right'); ax[1].set_title('(b) exact infinite armchair strip',fontsize=10)
plt.tight_layout(); plt.savefig('manuscript/figs/fig_honey_rule.pdf'); plt.savefig('fig_honey_rule.png',dpi=110)
