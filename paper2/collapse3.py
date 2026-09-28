import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
d=json.load(open('../data/paper2_ratio_fine.json')); base={r[0]:r for r in d['0']}
ts=[t for t in d if 0<float(t)<=0.1]
def curves(nu,which,Nmin=21):
    out=[]
    for t in ts:
        a=np.array([r for r in d[t] if r[0]>=Nmin]); N=a[:,0]
        r=a[:,which]/np.array([base[n][which] for n in N])
        out.append((t,np.log(N*float(t)**nu),np.log(r)))
    return out
def cost(nu,which):
    c=curves(nu,which); s=k=0
    for i,(_,xi,yi) in enumerate(c):
        for j,(_,xj,yj) in enumerate(c):
            if i==j: continue
            m=(xj>=xi.min())&(xj<=xi.max()); s+=np.sum((np.interp(xj[m],xi,yi)-yj[m])**2); k+=m.sum()
    return s/max(k,1)
nus=np.linspace(0.7,1.4,141)
fig,ax=plt.subplots(1,3,figsize=(13,3.8))
for which,name in [(1,'A (0,0)'),(2,'B (0,N-1)')]:
    c=np.array([cost(n,which) for n in nus]); b=nus[c.argmin()]; ok=nus[c<2*c.min()]
    print(name,'nu=%.3f  [%.3f,%.3f]  rms(best)=%.1e  rms(nu=1)=%.1e'%(b,ok.min(),ok.max(),np.sqrt(c.min()),np.sqrt(cost(1,which))))
    ax[0].semilogy(nus,np.sqrt(c),label=name)
    for t,x,y in curves(1.0,which,13):
        ax[which].plot(np.exp(x),np.exp(y),'o-',ms=3,label='t=%s'%t)
    xx=np.logspace(0.3,1,20); ax[which].plot(xx,0.83*xx**0.25,'k--',lw=1,label=r'$\propto (Nt)^{1/4}$')
    ax[which].set_xscale('log'); ax[which].set_yscale('log'); ax[which].set_xlabel('N t')
    ax[which].set_ylabel(r'$R(N,t)/R(N,0)$'); ax[which].set_title('corner '+name+r', $\nu=1$')
ax[0].set_xlabel(r'$\nu$'); ax[0].set_ylabel('collapse rms (log)'); ax[0].legend()
ax[1].legend(fontsize=7); plt.tight_layout(); plt.savefig('fig_collapse.png',dpi=110)
