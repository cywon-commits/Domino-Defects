import json,sys,numpy as np,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0,'.'); from symbol_xi import kappa
plt.rcParams.update({'font.size':14,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'#fbfaf6','axes.facecolor':'#fbfaf6','savefig.facecolor':'#fbfaf6'})
d=json.load(open('../data/paper2_pair_big.json')); cols=['#1f3a5f','#3d6a99','#6f9bc7','#e07b39','#c0501f','#8a3a12']
fig,ax=plt.subplots(1,2,figsize=(14,5.6))
for t,cl in zip(['0.1','0.15','0.2','0.3','0.4','0.5'],cols):
    C=np.array(d[t]['C']); r=C[:,0].astype(int); C=C[:,1]; k=kappa(float(t),2001)
    q=np.exp(-4*k); Ci=np.mean([C[-1]-(C[-3]-C[-1])*q/(1-q), C[-2]-(C[-4]-C[-2])*q/(1-q)])
    P=C-Ci; m=(r%2==1)&(np.abs(P)>1e-10*Ci); rr=r[m]; Y=np.abs(P[m])*np.exp(2*k*rr)
    s=np.diff(np.log(Y))/np.diff(np.log(rr)); rm=np.sqrt(rr[1:]*rr[:-1])
    ax[0].plot(rm,s,'o-',ms=4,color=cl,label='t = %s'%t)
    ax[1].semilogy(rr*2*k,Y*rr,'o-',ms=4,color=cl,label='t = %s'%t)
ax[0].axhline(-1,color='k',ls=':',lw=1.2); ax[0].set_ylim(-3.6,0); ax[0].set_xscale('log')
ax[0].set_xlabel('r (odd)'); ax[0].set_ylabel(r'$d\ln[(C-C_\infty)e^{2\kappa r}]/d\ln r$'); ax[0].legend(frameon=False,fontsize=11)
ax[0].set_title(r'local prefactor exponent $\to -1$ iff rate $=2\kappa$')
ax[1].set_xlabel(r'$2\kappa r$'); ax[1].set_ylabel(r'$(C-C_\infty)\,r\,e^{2\kappa r}$'); ax[1].set_title('amplitude A(t); vanishes at t = 1/2')
plt.tight_layout(); plt.savefig('fig_mono.png',dpi=140)
