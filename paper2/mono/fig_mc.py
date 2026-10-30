# Fig: interacting dimers MC. (a) pair weights, (b) edge vs bulk exponent, (c) corner/centre ratio vs L
import numpy as np, json, sys, os
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
S=sys.argv[1]; R=sys.argv[2]   # scratch, repo
W,H=64,128
summ=json.load(open(f'{R}/data/paper2_worm_summary.json'))
gcyl=dict(zip([float(j) for j in summ['J']],summ['bulk'])); ecyl=dict(zip([float(j) for j in summ['J']],summ['edge']))
def row(fn,r):
    h=np.fromfile(fn).reshape(H,W)[r]; x=np.arange(1,W,2); return (W/np.pi)*np.sin(np.pi*x/W),h[x]
cols=plt.cm.viridis(np.linspace(0,0.9,6))
fig,ax=plt.subplots(1,3,figsize=(13,3.9))
for i,J in enumerate(['1.0','0.5','0','-0.4','-0.7','-1.0']):
    lb,cb=row(f'{S}/mc/out/J{J}_bulk.bin',64); le,ce=row(f'{S}/mc/out/J{J}_edge.bin',0)
    m=lb<=20.4; ax[0].loglog(lb[m],cb[m]/cb[0],'o',color=cols[i],ms=4,label=r'$\beta v=%s$'%J)
    m=le<=20.4; ax[0].loglog(le[m],ce[m]/ce[0],'s',mfc='none',color=cols[i],ms=4)
ax[0].set_xlabel('chord distance'); ax[0].set_ylabel('pair weight (normalised)'); ax[0].legend(fontsize=7,ncol=2)
ax[0].set_title('(a) bulk (filled) and edge (open)',fontsize=10)
g=np.array(summ['bulk']); e=np.array(summ['edge'])
ax[1].plot(g,e,'o',color='#1f3a5f'); xx=np.linspace(0,2,10); ax[1].plot(xx,2*xx,'--',color='#c0501f',label='edge = 2 g')
ax[1].set_xlabel('bulk exponent g'); ax[1].set_ylabel('edge exponent'); ax[1].legend(); ax[1].set_title('(b) ratio 2.01 ± 0.05',fontsize=10)
def ratio(fn,L,r=2):
    h=np.fromfile(fn).reshape(L,L); x,y=np.meshgrid(range(L),range(L),indexing='ij'); b=((x+y)%2==0)
    c=n=0
    for cx,cy in [(0,0),(0,L-1),(L-1,0),(L-1,L-1)]:
        m=b&(abs(x-cx)+abs(y-cy)<=r); c+=h[m].sum(); n+=m.sum()
    mc=b&(abs(x-L//2)+abs(y-L//2)<=r); return (c/n)/(h[mc].sum()/mc.sum())
def ginterp(J):
    if J in gcyl: return gcyl[J]
    Js=sorted(gcyl); return float(np.interp(J,Js,[gcyl[j] for j in Js]))
out={}
for i,J in enumerate(['0.5','0','-0.5','-1.0']):
    Ls=[L for L in [17,25,33,49,65] if os.path.exists(f'{S}/mono/out/J{J}_L{L}.bin')]
    if len(Ls)<3: continue
    rr=np.array([ratio(f'{S}/mono/out/J{J}_L{L}.bin',L) for L in Ls]); Ls=np.array(Ls)
    p,cv=np.polyfit(np.log(Ls),np.log(rr),1,cov=True) if len(Ls)>3 else (np.polyfit(np.log(Ls),np.log(rr),1),np.zeros((2,2)))
    m=Ls>=25; p2=np.polyfit(np.log(Ls[m]),np.log(rr[m]),1)
    gg=ginterp(float(J)); out[J]=dict(L=Ls.tolist(),ratio=rr.tolist(),gamma_all=-p[0],err=float(np.sqrt(cv[0,0])),gamma_L25=-p2[0],g=gg,pred=gg/2)
    c=plt.cm.plasma(0.15+0.22*i)
    ax[2].loglog(Ls,rr,'o',color=c,label=r'$\beta v=%s$: $\gamma=%.2f$ (pred %.2f)'%(J,-p2[0],gg/2))
    ax[2].loglog(Ls,rr[-1]*(Ls/Ls[-1])**(-gg/2),'--',color=c,lw=1)
ax[2].set_xlabel('L (odd board)'); ax[2].set_ylabel('corner / centre'); ax[2].legend(fontsize=7); ax[2].set_title(r'(c) corner exponent vs $g/2$ (dashed)',fontsize=10)
plt.tight_layout(); plt.savefig(f'{R}/paper2/manuscript/figs/fig_mc.pdf'); plt.savefig(f'{S}/fig_mc.png',dpi=110)
json.dump(out,open(f'{R}/data/paper2_mono_corner.json','w'),indent=1)
for k,v in out.items(): print(k,v['L'],'gamma(all)=%.3f±%.3f gamma(L>=25)=%.3f  g/2=%.3f'%(v['gamma_all'],v['err'],v['gamma_L25'],v['pred']))
