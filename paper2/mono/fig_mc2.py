import numpy as np, json, sys
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter,FixedLocator,FixedFormatter
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
S=sys.argv[1]; R='../..'
W,H=64,128
summ=json.load(open(f'{R}/data/paper2_worm_summary.json')); mc=json.load(open(f'{R}/data/paper2_mono2.json'))
def row(fn,r):
    h=np.fromfile(fn).reshape(H,W)[r]; x=np.arange(1,W,2); return (W/np.pi)*np.sin(np.pi*x/W),h[x]
cols=['#1f3a5f','#3a6ea5','#2a9d8f','#e9a03b','#dd6b20','#8a3a12']; mks=['o','v','^','D','s','P']
fig,ax=plt.subplots(1,4,figsize=(14,3.4))
for i,J in enumerate(['1.0','0.5','0','-0.4','-0.7','-1.0']):
    lb,cb=row(f'{S}/mc/out/J{J}_bulk.bin',64); le,ce=row(f'{S}/mc/out/J{J}_edge.bin',0)
    m=lb<=20.4; ax[0].loglog(lb[m],cb[m]/cb[0],mks[i],color=cols[i],ms=3.5,label=r'$\beta v=%s$'%J)
    m=le<=20.4; ax[0].loglog(le[m],ce[m]/ce[0],mks[i],mfc='none',color=cols[i],ms=3.5)
ax[0].set_xlabel('chord distance'); ax[0].set_ylabel('pair weight'); ax[0].legend(frameon=False,fontsize=7,ncol=2,loc='lower left')
ax[0].set_title('(a) bulk (filled), edge (open)',fontsize=9)
g=np.array(summ['bulk']); e=np.array(summ['edge'])
ax[1].plot(g,e,'o',color=cols[0]); xx=np.linspace(0,2,10); ax[1].plot(xx,2*xx,'--',color=cols[4],lw=1,label='edge = 2 × bulk')
ax[1].set_xlabel('bulk exponent g'); ax[1].set_ylabel('edge exponent'); ax[1].legend(frameon=False,fontsize=8); ax[1].set_title('(b) pair exponents',fontsize=9)
gc=dict(zip([float(j) for j in summ['J']],summ['bulk']))
def gi(J): Js=sorted(gc); return float(np.interp(J,Js,[gc[j] for j in Js]))
for i,J in enumerate(['0.5','0','-0.5','-1.0']):
    d=mc['corner'].get(J); 
    if not d: continue
    L=np.array(d['L']); r=np.array(d['ratio']); er=np.array(d['err']); gg=gi(float(J)); c=cols[[0,2,4,5][i]]
    ax[2].errorbar(L,r,yerr=er,fmt=mks[i],color=c,ms=4,capsize=2,label=r'$\beta v=%s$: $\gamma=%.2f\pm%.2f$ ($g/2=%.2f$)'%(J,d['gamma'],d['gerr'],gg/2))
    ax[2].plot(L,r[-1]*(L/L[-1])**(-gg/2),'--',color=c,lw=0.8)
ax[2].set_xscale('log'); ax[2].set_yscale('log'); ax[2].set_xlabel('L (odd board)'); ax[2].set_ylabel('corner / centre'); ax[2].legend(frameon=False,fontsize=6.5,loc='lower left'); ax[2].set_title('(c) corner exponent',fontsize=9)
for a in ax[2:]:
    a.xaxis.set_minor_formatter(NullFormatter()); a.xaxis.set_major_locator(FixedLocator([16,32,64])); a.xaxis.set_major_formatter(FixedFormatter(['16','32','64']))
for k,(key,lab,c,mk) in enumerate([('edge_-1.0','edge, g = 1.88',cols[5],'s'),('bulk_-1.0','interior, g = 1.88',cols[5],'o'),('edge_-0.4','edge, g = 0.84',cols[2],'s'),('bulk_-0.4','interior, g = 0.84',cols[2],'o')]):
    d=mc['vac'].get(key)
    if not d: continue
    L=np.array(d['L']); P=np.array(d['P']); e=np.array(d['err']); Rr=(1-P)/P; Re=e/P**2
    ax[3].errorbar(L,Rr,yerr=Re,fmt=mk,color=c,mfc=c if mk=='o' else 'none',ms=4,capsize=2,label=lab,ls='-',lw=0.8)
ax[3].set_xscale('log'); ax[3].set_yscale('log'); ax[3].set_xlabel('L'); ax[3].set_ylabel(r'$R=(1-P)/P$'); ax[3].legend(frameon=False,fontsize=7,loc='upper left'); ax[3].set_title('(d) vacancy: far/near weight',fontsize=9)
for a in ax[2:]:
    a.xaxis.set_minor_formatter(NullFormatter()); a.xaxis.set_major_locator(FixedLocator([16,32,64])); a.xaxis.set_major_formatter(FixedFormatter(['16','32','64']))
plt.tight_layout(); plt.savefig(f'{R}/paper2/manuscript/figs/fig_mc.pdf'); plt.savefig(f'{S}/fig_mc2.png',dpi=110)
