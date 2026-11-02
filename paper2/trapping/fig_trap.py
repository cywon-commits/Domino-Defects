import json,numpy as np,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter,FixedLocator,FixedFormatter
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
d=json.load(open('../../data/paper2_trapping.json')); C1,C2='#1f3a5f','#dd6b20'
cases=[('bulk_q1','interior, q = 1',C1,'o',-0.5),('edge_q1','edge, q = 1',C1,'s',-0.75),('bulk_q3','interior, q = 3',C2,'o',-1.5),('edge_q3','edge, q = 3',C2,'s',-2.75)]
fig,ax=plt.subplots(1,2,figsize=(7.2,3.0))
for k,lab,c,mk,g in cases:
    Ls=np.array(sorted(int(x) for x in d[k])); P=np.array([d[k][str(L)]['P_near3'] for L in Ls])
    ax[0].loglog(Ls,P,mk,color=c,mfc=c if 'bulk' in k else 'none',ms=5,label=r'%s  ($\gamma=%s$)'%(lab,{-0.5:'-1/2',-0.75:'-3/4',-1.5:'-3/2',-2.75:'-11/4'}[g]))
    if g>-2: ax[0].loglog(Ls,P[-1]*(Ls/Ls[-1])**(-(2+g)),'-',color=c,lw=0.8)
    pr=np.array(d[k]['80']['prof']); r,rho=pr[0],pr[1]; m=r>0
    ax[1].loglog(r[m],rho[m],mk,color=c,mfc=c if 'bulk' in k else 'none',ms=3,alpha=0.8)
ax[0].axhline(0.34,color=C2,ls='--',lw=0.8); ax[0].text(17,0.36,r'$P_\infty\approx0.34$',color=C2,fontsize=8)
ax[0].set_xlabel('L'); ax[0].set_ylabel('P (defect within distance 3)'); ax[0].legend(frameon=False,fontsize=7,loc='lower left')
ax[0].xaxis.set_minor_formatter(NullFormatter()); ax[0].xaxis.set_major_locator(FixedLocator([16,32,64])); ax[0].xaxis.set_major_formatter(FixedFormatter(['16','32','64']))
ax[0].set_title(r'(a) lines: $L^{-(2+\gamma)}$',fontsize=9)
ax[1].set_xlabel('distance d from the fixed charge'); ax[1].set_ylabel(r'$\rho(d)$, L = 80'); ax[1].set_title('(b) profiles (same symbols)',fontsize=9)
plt.tight_layout(); plt.savefig('../manuscript/figs/fig_trap.pdf'); plt.savefig('fig_trap.png',dpi=110)
