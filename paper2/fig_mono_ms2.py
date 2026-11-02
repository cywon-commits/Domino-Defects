import json,numpy as np,matplotlib,matplotlib.ticker; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':8.5,'axes.spines.top':False,'axes.spines.right':False})
d=json.load(open('../data/paper2_pair_big.json')); cols=['#1f3a5f','#3d6a99','#6f9bc7','#e07b39','#c0501f','#8a3a12']
def rate(t):
    if t==0.5: return np.log(2+np.sqrt(3))
    mu=np.sqrt(1-4*t*t); return -np.log(2-mu-2*np.sqrt(1-t*t-mu))
Cinf=lambda t:0.5*np.sqrt(t/(2*t*(2+t*t)+(1+2*t*t)*np.sqrt(2+t*t)))
fig,ax=plt.subplots(1,2,figsize=(5.8,2.6))
for t,cl in zip(['0.1','0.15','0.2','0.3','0.4','0.5'],cols):
    tt=float(t); C=np.array(d[t]['C']); r=C[:,0].astype(int); C=C[:,1]; P=C-Cinf(tt); k=rate(tt)
    m=(r%2==1)&(np.abs(P)>2e-10); rr=r[m]; Y=np.abs(P[m])*np.exp(k*rr)
    s=np.diff(np.log(Y))/np.diff(np.log(rr)); rm=np.sqrt(rr[1:]*rr[:-1])
    ax[0].plot(rm,s,'o-',ms=4,color=cl,label='t = %s'%t); ax[1].semilogy(rr*k,Y*rr,'o-',ms=4,color=cl,label='t = %s'%t)
ax[0].axhline(-1,color='k',ls=':',lw=1.2); ax[0].axhline(0,color='k',ls=':',lw=0.8); ax[0].set_ylim(-1.8,0.3); ax[0].set_xscale("log"); ax[0].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax[0].set_xlabel('r (odd)'); ax[0].set_ylabel(r'$d\ln[(C-C_\infty)e^{2r\ln\eta_2}]/d\ln r$'); ax[0].legend(frameon=False,fontsize=6.5)
ax[0].set_title(r'(b) prefactor exponent: $-1$ for $t<1/2$, $0$ at $t=1/2$',fontsize=7.5)
ax[1].set_xlabel(r'$2r\ln\eta_2$'); ax[1].set_ylabel(r'$(C-C_\infty)\,r\,e^{2r\ln\eta_2}$'); ax[1].set_title(r'(c) $r^{-1}$ amplitude; grows linearly at $t=1/2$',fontsize=7.5)
plt.tight_layout(); plt.savefig('manuscript/figs/fig03b_mono.pdf'); plt.savefig('/tmp/f3b.png',dpi=60)
