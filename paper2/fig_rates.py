# Sec. VI figure: (a) decay rates of C-C_inf per unit length vs t; (b) a2(N)/N^2
import json,numpy as np,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
D='../data/'; C1,C2,C3='#1f3a5f','#dd6b20','#8a8f96'
def eta2(t):
    mu=np.sqrt(1-4*t*t+0j); x=2-mu-2*np.sqrt(1-t*t-mu); return np.real(-0.5*np.log(x))
fig,ax=plt.subplots(1,2,figsize=(10,3.9))
tt=np.linspace(0.02,0.5,200); ax[0].plot(tt,2*eta2(tt),'-',color=C1,lw=2,label=r'$2\ln\eta_2$ (Basor–Bleher), row')
ax[0].plot(tt,2*tt,':',color=C3,lw=1.2,label=r'$2t$')
o=json.load(open(D+'paper2_antidiag.json')); o96=json.load(open(D+'paper2_antidiag96.json')); out={}
for t in ['0.2','0.3','0.4','0.5']:
    d=np.array(o96[t]) if t in o96 else np.array(o[t]['res']['anti']['d'][:10])
    r=np.log(d[6]/d[8])/2/np.sqrt(2)   # odd k=7,9, per Euclidean length
    out[t]=r; ax[0].plot(float(t),r,'s',color=C2,ms=7,mfc='none',mew=1.6,label='anti-diagonal (1,−1), measured' if t=='0.2' else None)
ax[0].plot(0.5,np.log(2+np.sqrt(3)),'o',color=C1,ms=7,label=r'row and (1,1) at $t=1/2$: $\ln(2+\sqrt{3})$')
ax[0].set_xlabel('t'); ax[0].set_ylabel('decay rate per unit length'); ax[0].legend(frameon=False,fontsize=8.5,loc='upper left'); ax[0].set_title('(a) screening rates',fontsize=11)
a=np.array(json.load(open(D+'paper2_deriv.json'))); N=a[:,0]; m=N>=29; NN=np.linspace(25,110,200)
meas101={2:203.449,4:208.065}
for k_,n,cl in [(2,'A',C1),(4,'B',C2)]:
    X=lambda x:np.c_[x**2,x*np.log(x),np.ones_like(x)]; c,*_=np.linalg.lstsq(X(N[m]),a[m,k_],rcond=None)
    ax[1].plot(N,a[:,k_]/N**2,'o',ms=6,color=cl,label='corner %s'%n); ax[1].plot(NN,(X(NN)@c)/NN**2,'-',lw=1.5,color=cl)
    ax[1].plot([101],[meas101[k_]/101**2],'*',ms=14,color=cl,mec='k')
ax[1].set_xlabel('N'); ax[1].set_ylabel(r'$a_2(N)/N^2$'); ax[1].legend(frameon=False,fontsize=9,loc='upper right')
ax[1].set_title(r'(b) fit $cN^2+eN\ln N+f$ on $29\leq N\leq81$; stars: $N=101$',fontsize=10)
plt.tight_layout(); plt.savefig('manuscript/figs/fig_rates.pdf'); plt.savefig('fig_rates.png',dpi=110)
print({k:round(v,3) for k,v in out.items()}, {t:round(2*eta2(float(t)),3) for t in out})
