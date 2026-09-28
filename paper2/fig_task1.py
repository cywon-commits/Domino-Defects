import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, sys
sys.path.insert(0,'.'); from symbol_xi import kappa
f=json.load(open('../data/paper2_ratio_fine.json')); p=json.load(open('../data/paper2_pair_bulk.json'))
fig,ax=plt.subplots(1,3,figsize=(13.5,3.9))
for t in ['0','0.01','0.02','0.05','0.1','0.2']:
    a=np.array(f[t]); ax[0].plot(a[:,0],a[:,1],'o-',ms=3,label='t=%s'%t)
ax[0].set_xscale('log'); ax[0].set_yscale('log'); ax[0].set_xlabel('N'); ax[0].set_ylabel('R = P(corner A)/P(center)')
ax[0].set_title('(a) forced single defect, odd N×N'); ax[0].legend(fontsize=7)
for t in ['0','0.1','0.2','0.5','1.0']:
    a=np.array(p[t]); ax[1].plot(a[:,0],a[:,1],'o-',ms=3,label='t=%s'%t)
ax[1].set_xscale('log'); ax[1].set_yscale('log'); ax[1].set_xlabel('r'); ax[1].set_ylabel('C(r)')
ax[1].set_title('(b) bulk pair, 64×64: deconfinement'); ax[1].legend(fontsize=7)
tt=np.linspace(0.01,1,60); k=np.array([kappa(x,801) for x in tt])
ax[2].plot(tt,k*1,'-',label=r'exact $\kappa(t)$ from symbol'); ax[2].plot(tt,np.arctanh(np.minimum(tt,0.999)),'--',label=r'artanh$(t)$')
ax[2].set_ylim(0,1.2); ax[2].set_xlabel('t'); ax[2].set_ylabel(r'$1/\xi_f$ (per site)'); ax[2].set_title(r'(c) fermion gap: $\xi_f=1/\mathrm{artanh}\,t$'); ax[2].legend(fontsize=8)
plt.tight_layout(); plt.savefig('fig_task1_summary.png',dpi=120)
