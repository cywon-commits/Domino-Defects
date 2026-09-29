import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, sys
sys.path.insert(0,'.'); from symbol_xi import kappa
plt.rcParams.update({'font.size':15,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':1.2,
    'figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white'})
cols=['#1f3a5f','#3d6a99','#6f9bc7','#e07b39','#c0501f','#8a3a12']
D='../data/'
# F1 lattice geometry
fig,ax=plt.subplots(figsize=(6,6)); N=7
for r in range(N):
    for c in range(N):
        if c+1<N: ax.plot([c,c+1],[r,r],color='#1f3a5f',lw=2)
        if r+1<N: ax.plot([c,c],[r,r+1],color='#1f3a5f',lw=2)
        if r+1<N and c+1<N: ax.plot([c,c+1],[r,r+1],color='#e07b39',lw=2,ls='--')
xs,ys=np.meshgrid(range(N),range(N)); ax.scatter(xs,ys,s=60,color='#1f3a5f',zorder=3)
ax.scatter([3],[3],s=380,facecolor='#fbfaf6',edgecolor='#c0501f',lw=3,zorder=4)
ax.plot([3.12,3.12],[3,-0.6],color='#c0501f',lw=2.5,ls=':')
ax.text(-0.35,-0.45,'A',fontsize=20,color='#1f3a5f',weight='bold'); ax.text(N-0.8,-0.45,'B',fontsize=20,color='#1f3a5f',weight='bold')
ax.text(3.3,1.4,'string',color='#c0501f',fontsize=15)
ax.set_aspect('equal'); ax.axis('off'); plt.tight_layout(); plt.savefig('manuscript/figs/fig01a_lattice.pdf',dpi=200)
# F2 corner ratio
f=json.load(open(D+'paper2_ratio_fine.json'))
fig,ax=plt.subplots(figsize=(8.5,6))
for t,cl in zip(['0','0.01','0.02','0.05','0.1','0.2'],cols):
    a=np.array(f[t]); ax.plot(a[:,0],a[:,1],'o-',ms=6,lw=2,color=cl,label='t = %s'%t)
NN=np.array([13,81]); ax.plot(NN,1.84*NN**-0.25,'k:',lw=1.5,label=r'$\propto N^{-1/4}$')
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('N  (odd N x N board)'); ax.set_ylabel('R = P(corner A) / P(center)')
ax.legend(frameon=False,fontsize=13,ncol=2); plt.tight_layout(); plt.savefig('manuscript/figs/fig02_corner_ratio.pdf',dpi=200)
# F3 pair correlation
p=json.load(open(D+'paper2_pair_bulk.json'))
fig,ax=plt.subplots(figsize=(8.5,6))
for t,cl in zip(['0','0.1','0.2','0.5','1.0'],[cols[0],cols[2],cols[3],cols[4],cols[5]]):
    a=np.array(p[t]); m=a[:,1]>0; ax.plot(a[m,0],a[m,1],'o-',ms=5,lw=2,color=cl,label='t = %s'%t)
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('separation r'); ax.set_ylabel('C(r) = Z(two defects) / Z')
ax.legend(frameon=False,fontsize=13); plt.tight_layout(); plt.savefig('manuscript/figs/fig03a_pair.pdf',dpi=200)
# F4 kappa
tt=np.linspace(0.005,1,120); k=np.array([kappa(x,1201) for x in tt])
fig,ax=plt.subplots(figsize=(8.5,6))
ax.plot(tt,k,'-',lw=3,color=cols[0],label=r'exact $\kappa(t)$ from Kasteleyn symbol')
ax.plot(tt,np.arctanh(np.minimum(tt,0.999)),'--',lw=2.5,color=cols[3],label=r'artanh$(t)$')
ax.set_ylim(0,1.2); ax.set_xlabel('diagonal weight t'); ax.set_ylabel(r'$1/\xi_f$  (per lattice site)')
ax.legend(frameon=False,fontsize=13); plt.tight_layout(); plt.savefig('manuscript/figs/fig04a_kappa.pdf',dpi=200)
# F5 a2
a=np.array(json.load(open(D+'paper2_deriv.json'))); N=a[:,0]; m=N>=29
fig,ax=plt.subplots(figsize=(8.5,6)); NN=np.linspace(25,110,200)
meas101={'A':203.449,'B':208.065}
for k_,n,cl in [(2,'A',cols[0]),(4,'B',cols[3])]:
    X=lambda x:np.c_[x**2,x*np.log(x),np.ones_like(x)]; c,*_=np.linalg.lstsq(X(N[m]),a[m,k_],rcond=None)
    ax.plot(N,a[:,k_]/N**2,'o',ms=7,color=cl,label='corner %s'%n); ax.plot(NN,(X(NN)@c)/NN**2,'-',lw=1.8,color=cl)
    ax.plot([101],[meas101[n]/101**2],'*',ms=18,color=cl,mec='k')
ax.set_xlabel('N'); ax.set_ylabel(r'$a_2(N)/N^2$'); ax.legend(frameon=False,fontsize=13)
ax.text(62,0.0395,'lines: fit  $cN^2+eN\\ln N+f$  on N ≤ 81\nstars: N = 101 (not used in fit)',fontsize=13)
plt.tight_layout(); plt.savefig('manuscript/figs/fig04b_a2.pdf',dpi=200)
# F6 t=1/2
from pftri import logZ, logZ_multi
L=48; t=0.5; lz=logZ(L,L,t); c0=L//2; q=2-np.sqrt(3)
fig,ax=plt.subplots(figsize=(8.5,6))
for name,step,cl,mk in [('(1,0)',(1,0),cols[0],'o'),('(0,1)',(0,1),cols[2],'s'),('(1,1)',(1,1),cols[3],'^'),('(1,2)',(1,2),cols[5],'D'),('(1,-1)',(1,-1),'#888888','x')]:
    ds=[];vals=[]
    for kk in range(1,9 if step!=(1,2) else 5):
        dr,dc=step[0]*kk,step[1]*kk; u=(c0-dr//2,c0-dc//2); v=(u[0]+dr,u[1]+dc)
        C=np.exp(logZ_multi(L,L,t,[u,v])-lz); d=abs(dr)+abs(dc)
        if abs(C-1/6)>1e-12: ds.append(d); vals.append(abs(C-1/6))
    ax.semilogy(ds,vals,mk,ms=9,color=cl,label='direction '+name,mfc='none' if mk!='x' else None,mew=2)
dd=np.arange(1,17); ax.semilogy(dd,q**dd/6,'k-',lw=1.2,label=r'$(2-\sqrt{3})^{d}/6$')
ax.set_xlabel(r'$d=|\Delta r|+|\Delta c|$'); ax.set_ylabel(r'$|C-1/6|$'); ax.legend(frameon=False,fontsize=12)
plt.tight_layout(); plt.savefig('manuscript/figs/fig05_half.pdf',dpi=200)
