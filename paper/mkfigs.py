import json, math, pickle, numpy as np, os, sys
import sys; sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','lattice'))
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon, Circle, FancyArrowPatch
from matplotlib.colors import LinearSegmentedColormap
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','data'))
plt.rcParams.update({'font.family':'serif','font.size':8,'axes.labelsize':8,'axes.titlesize':8,'legend.fontsize':7,'xtick.labelsize':7,'ytick.labelsize':7,'mathtext.fontset':'cm','axes.linewidth':0.6,'lines.linewidth':1.1})
C1,C2,C3,C4='#1f5f8b','#c0562b','#3a8f5c','#7a5195'
D=json.load(open('figdata.json'))
cmap=LinearSegmentedColormap.from_list('d',['#f7f7f7','#9ecae1','#2171b5','#08306b'])
def save(fig,name): fig.savefig(f'../paper/figs/{name}.pdf',bbox_inches='tight'); plt.close(fig); print('saved',name)
def heat(ax,M,vmax,title,vmin=0):
    im=ax.imshow(np.array(M),origin='lower',cmap=cmap,vmin=vmin,vmax=vmax,interpolation='nearest'); ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title); return im
def mark_removed(ax,x,y,label,dx,dy,ha):
    ax.plot(x,y,marker='s',ms=5.5,mfc=C2,mec='k',mew=0.5,clip_on=False,zorder=6)
    ax.text(x+dx,y+dy,label,color=C2,fontsize=6.5,ha=ha,va='center',clip_on=False)
def loglog_series(ax,xs,ys,color,label,slope=None,anchor=-1,marker='o'):
    ax.loglog(xs,ys,marker,color=color,ms=3.5,label=label)
    if slope is not None:
        xa=np.array([min(xs)*0.9,max(xs)*1.1]); y0=ys[anchor]*(xa/xs[anchor])**slope; ax.loglog(xa,y0,'--',color=color,lw=0.8)
# ---------- Fig 1 schematic ----------
fig,axs=plt.subplots(1,2,figsize=(7,3.1))
ax=axs[0]; N=8
for y in range(N):
    for x in range(N):
        col='#333333' if (x+y)%2==0 else '#e8e8e8'
        if (x,y)==(0,0): ax.add_patch(Rectangle((x,y),1,1,fc='white',ec='#c0562b',hatch='////',lw=0.8)); continue
        ax.add_patch(Rectangle((x,y),1,1,fc=col,ec='#999999',lw=0.3))
ax.add_patch(Circle((4.5,5.5),0.28,fc=C2,ec='k',lw=0.5)); ax.text(4.95,5.75,r'defect $e=-1$',fontsize=7)
for (px,py,txt,col) in [(-0.35,-0.35,r'$+\frac{3}{4}$',C2),(8.1,-0.35,r'$+\frac{1}{4}$',C1),(-0.55,8.1,r'$+\frac{1}{4}$',C1),(8.1,8.1,r'$-\frac{1}{4}$',C1)]:
    ax.text(px,py,txt,color=col,fontsize=9,ha='center')
ax.set_xlim(-1,9.2); ax.set_ylim(-1,9.2); ax.set_aspect('equal'); ax.axis('off'); ax.set_title('(a) corner charges $q_c=-J_c/4$')
ax=axs[1]
ax.plot([0,6],[0,0],'k',lw=1.2); ax.plot([0,0],[0,4.2],'k',lw=1.2)
ax.fill_between([0,6],0,-0.25,color='#dddddd'); ax.fill_betweenx([0,4.2],0,-0.25,color='#dddddd')
z=(3.2,2.2)
for (x,y,lab,al) in [(z[0],z[1],'defect',1.0),(z[0],-z[1],'image',0.45),(-z[0],z[1],'image',0.45),(-z[0],-z[1],'image',0.3)]:
    ax.add_patch(Circle((x,y),0.22,fc=C2,ec='k',lw=0.5,alpha=al)); ax.text(x+0.3,y+0.15,lab,fontsize=7,alpha=max(al,0.6))
ax.add_patch(Circle((0,0),0.18,fc=C1,ec='k',lw=0.5)); ax.text(0.25,0.25,r'corner $q_c$',fontsize=7,color=C1)
ax.annotate('',xy=(0.5,0.4),xytext=(z[0]-0.25,z[1]-0.2),arrowprops=dict(arrowstyle='->',color=C1,lw=0.8))
ax.text(1.2,1.55,'Coulomb\nattraction',fontsize=6.5,color=C1,rotation=33)
ax.text(3.4,3.5,r'$\rho\propto|f^\prime|^{1/4}(\mathrm{Im}f)^{1/4}\prod_c|f-f_c|^{e q_c}$',fontsize=7)
ax.set_xlim(-4,6.3); ax.set_ylim(-3,4.4); ax.set_aspect('equal'); ax.axis('off'); ax.set_title('(b) same-sign images at the edge')
save(fig,'fig1_schematic')
# ---------- Fig 2 single corner ----------
fig=plt.figure(figsize=(7,2.35)); gs=fig.add_gridspec(1,4,width_ratios=[1,1,0.06,1.5],wspace=0.12)
a0=fig.add_subplot(gs[0]); a1=fig.add_subplot(gs[1]); cax=fig.add_subplot(gs[2]); a2=fig.add_subplot(gs[3])
im=heat(a0,D['single64']['data'],2.0,'(a) exact, $N=64$'); heat(a1,D['single64']['theory'],2.0,'(b) Coulomb gas')
for _a in (a0,a1): mark_removed(_a,-0.5,-0.5,'a1',0.9,-1.6,'left')
cb=fig.colorbar(im,cax=cax); cb.set_label(r'$\rho/\rho_{\rm unif}$')
a2.set_position([a2.get_position().x0+0.15,a2.get_position().y0,a2.get_position().width-0.15,a2.get_position().height])
Ns=sorted(int(k) for k in D['corner']); arr=np.array([D['corner'][str(n)] for n in Ns])
loglog_series(a2,Ns,arr[:,0],C2,r'cell next to removed a1 (slope $3/4$)',0.75)
loglog_series(a2,Ns,arr[:,1],C1,r'corner h1 (slope $-1/4$)',-0.25)
loglog_series(a2,Ns,arr[:,2],C3,r'cell next to h8 (slope $-5/4$)',-1.25)
a2.set_xlabel('board size $N$'); a2.set_ylabel(r'$\rho/\rho_{\rm unif}$'); a2.set_ylim(5e-3,4e2); a2.legend(frameon=False,loc='upper left',fontsize=6); a2.set_title('(c) corner scaling (dashed: theory)')
save(fig,'fig2_single')
# ---------- Fig 3 mutilated ----------
fig=plt.figure(figsize=(7,4.3)); gs=fig.add_gridspec(2,4,width_ratios=[1,1,0.05,1.25],hspace=0.5,wspace=0.18)
a=fig.add_subplot(gs[0,0]); b=fig.add_subplot(gs[0,1]); cax=fig.add_subplot(gs[0,2]); c=fig.add_subplot(gs[0,3])
im=heat(a,D['mut96f']['data'],2.0,'(a) exact, $N=96$'); heat(b,D['mut96f']['theory'],2.0,'(b) Coulomb gas')
nb=len(D['mut96f']['data'])
for _a in (a,b): mark_removed(_a,-0.5,-0.5,'a1',1.5,-2.6,'left'); mark_removed(_a,nb-0.5,nb-0.5,'h8',-1.5,2.6,'right')
fig.colorbar(im,cax=cax).set_label(r'$\rho/\rho_{\rm unif}$')
c.set_position([c.get_position().x0+0.14,c.get_position().y0,c.get_position().width-0.14,c.get_position().height])
bins=np.array(D['mut96']['bins']); mid=(bins[1:]+bins[:-1])/2; w=bins[1]-bins[0]
c.bar(mid,D['mut96']['hd'],width=w*0.85,color=C1,alpha=0.55,label='exact'); c.plot(mid,D['mut96']['ht'],'o',color=C2,ms=3.5,label='theory')
c.set_xlabel('defect separation $/N$'); c.set_ylabel('probability'); c.set_ylim(0,0.19); c.legend(frameon=False,loc='upper right',fontsize=6,handlelength=1.2); c.set_title('(c) two-defect separation')
d=fig.add_subplot(gs[1,0:3]); mc=np.array(D['mutcorner'])
loglog_series(d,mc[:,0],mc[:,1],C2,r'next to removed a1 ($+3/4$)',0.75); loglog_series(d,mc[:,0],mc[:,2],C1,r'corner h1 ($-1/4$)',-0.25)
d.set_xlabel('$N$'); d.set_ylabel(r'$\rho/\rho_{\rm unif}$'); d.legend(frameon=False); d.set_title('(d) corner scaling, $N=8$–$96$')
e=fig.add_subplot(gs[1,3]); mp_=D['mutprof']; e.plot(mp_['r'],mp_['d'],'o',color=C1,ms=3.5,label='exact'); e.plot(mp_['r'],mp_['t'],'-',color=C2,label='theory')
e.set_xscale('log'); e.set_xlabel('distance $r$'); e.set_ylabel('conditional density'); e.legend(frameon=False); e.set_title('(e) around a central defect')
save(fig,'fig3_mutilated')
# ---------- Fig 4 odd board ----------
O=D['odd63s']; fig=plt.figure(figsize=(7,4.3)); gs=fig.add_gridspec(2,4,width_ratios=[1,1,1,0.05],hspace=0.45,wspace=0.12)
from matplotlib.ticker import NullFormatter, FixedLocator
a=fig.add_subplot(gs[0,0]); b=fig.add_subplot(gs[0,1]); c=fig.add_subplot(gs[0,2]); cax=fig.add_subplot(gs[0,3])
lo,hi=0.5,1.3
im=heat(a,np.array(O['ee'],dtype=float),hi,'(a) exact, (even,even) cells',lo)
heat(b,np.array(O['oo'],dtype=float),hi,'(b) exact, (odd,odd) cells',lo)
heat(c,np.array(O['tee'],dtype=float),hi,r'(c) $\rho\propto R^{1/4}$',lo)
fig.colorbar(im,cax=cax).set_label(r'$\rho/\rho_{\rm unif}$')
d=fig.add_subplot(gs[1,0:2]); cut=np.array(O['cut']); ev=cut[cut[:,0]%2==0]; od=cut[cut[:,0]%2==1]
d.plot(ev[:,0],ev[:,1],'o',ms=3,color=C1,label='exact, (even,even)'); d.plot(od[:,0],od[:,1],'s',ms=3,color=C2,label='exact, (odd,odd)')
d.plot(cut[:,0],cut[:,2],'-',color='k',lw=0.9,label=r'$\rho\propto R^{1/4}$')
d.set_xlim(-1,63); d.set_xlabel('row $y$ along the central column'); d.set_ylabel(r'$\rho/\rho_{\rm unif}$'); d.legend(frameon=False,loc='lower center'); d.set_title('(d) cut through the centre, $N=63$')
e=fig.add_subplot(gs[1,2:4]); sq=np.array(D['triratio']['sq']); loglog_series(e,sq[:,0],sq[:,1],C3,'corner/centre',-0.25)
e.set_position([e.get_position().x0+0.07,e.get_position().y0,e.get_position().width-0.07,e.get_position().height])
e.yaxis.set_major_locator(FixedLocator([0.7,0.8,0.9,1.0,1.1])); e.yaxis.set_minor_formatter(NullFormatter()); e.set_yticklabels(['0.7','0.8','0.9','1.0','1.1'])
e.xaxis.set_major_locator(FixedLocator([7,15,31,63])); e.xaxis.set_minor_formatter(NullFormatter()); e.set_xticklabels(['7','15','31','63'])
e.set_xlabel('$N$ (odd)'); e.set_ylabel('corner/centre ratio'); e.legend(frameon=False); e.set_title('(e) slope $-1/4$')
save(fig,'fig4_odd')
# ---------- Fig 5 pairs, symbol, amplitudes ----------
from mpmath import barnesg
G=float(barnesg(0.5))
fig,axs=plt.subplots(1,3,figsize=(7,2.3),gridspec_kw=dict(wspace=0.38))
bd=np.array(D['bnd'])[:6]
loglog_series(axs[0],bd[:,0],bd[:,1],C2,'boundary pair ($-1$)',-1,3); loglog_series(axs[0],bd[:,0],bd[:,2],C1,'bulk pair ($-1/2$)',-0.5,3)
axs[0].set_xlabel('separation $d$'); axs[0].set_ylabel('pair weight'); axs[0].legend(frameon=False,loc='lower left'); axs[0].set_title('(a) monomer pairs, $160\\times160$')
ph=np.linspace(-np.pi,np.pi,2001); c_=np.cos(ph); sig=np.sign(c_)*(c_+1j)/np.sqrt(1+c_**2)
axs[1].plot(ph/np.pi,np.angle(sig)/np.pi,'.',ms=0.8,color=C4); axs[1].set_xlabel(r'$\varphi/\pi$'); axs[1].set_ylabel(r'$\arg\sigma/\pi$')
axs[1].set_title(r'(b) scalar symbol $\sigma(\varphi)$'); axs[1].axvline(0.5,ls=':',c='gray'); axs[1].axvline(-0.5,ls=':',c='gray')
def load(f): d=json.load(open(f)); return np.array(sorted((int(k),float(v)) for k,v in d.items()))
b=load('exact_bulk_15_1023.json'); h=load('exact_half_15_1023.json'); q=load('exact_quad_a.json'); qb=load('exact_quad.json'); q=np.vstack([q,qb])
Cb=math.pi*2**-0.75*G**4; Ch=2**0.375*G**2; Cq=2**(9/8)*G**2
axs[2].plot(1/b[:,0],b[:,1]*b[:,0]**0.5/Cb,'o-',ms=3,color=C1,label='bulk pair')
axs[2].plot(1/h[:,0],h[:,1]*h[:,0]**0.75/Ch,'s-',ms=3,color=C2,label='half plane')
rq=q[:,0]; rho=np.sqrt(rq**2+(rq+1)**2)
axs[2].plot(1/rq,q[:,1]*rho**0.75/Cq,'^-',ms=3,color=C3,label='quadrant')
axs[2].axhline(1,c='k',lw=0.6); axs[2].set_xscale('log'); axs[2].set_xlabel('$1/r$'); axs[2].set_ylabel('scaled weight / closed form'); axs[2].legend(frameon=False); axs[2].set_title('(c) amplitudes $\\to$ Barnes-$G$ forms')
save(fig,'fig5_pairs')
# ---------- Fig 6 Aztec ----------
fig,axs=plt.subplots(1,3,figsize=(7,2.35),gridspec_kw=dict(width_ratios=[1,1.2,1.2],wspace=0.35))
az=D['az_single']; n=az['n']; ax=axs[0]
from kreg import aztec as _az
from matplotlib.colors import Normalize
cmz=plt.get_cmap('YlGnBu'); nz=Normalize(-6,2)
for r,cc,lv in az['log10']:   # each defect-sublattice cell drawn as a diamond reaching its four neighbours (no checkerboard gaps)
    x,y=cc+0.5,r+0.5; ax.add_patch(Polygon([(x-1,y),(x,y+1),(x+1,y),(x,y-1)],closed=True,fc=cmz(nz(lv)),ec='none'))
_reg=set(tuple(map(int,v)) for v in _az(n))
for r in range(-2,2*n+2):
    for cc in range(-2,2*n+2):
        if (r,cc) not in _reg: ax.add_patch(Rectangle((cc,r),1,1,fc='white',ec='none',zorder=3))
ax.set_xlim(-7,2*n); ax.set_ylim(-1,2*n+1); ax.set_aspect('equal'); ax.axis('off'); ax.set_title(f'(a) one defect, $n={n}$')
th=np.linspace(0,2*np.pi,200); ax.plot(n+n/np.sqrt(2)*np.cos(th),n+n/np.sqrt(2)*np.sin(th),'--',c=C2,lw=0.7,zorder=4)
ax.add_patch(Rectangle((0,n),1,1,fc=C2,ec='k',lw=0.4,zorder=6)); ax.add_patch(Circle((0.5,n+0.5),2.2,fc='none',ec=C2,lw=0.9,zorder=6))
ax.annotate('removed\ntip',xy=(-1.6,n+1.6),xytext=(-6.5,n+8),color=C2,fontsize=6,ha='center',arrowprops=dict(arrowstyle='-',color=C2,lw=0.6))
cax_=ax.inset_axes([0.12,-0.10,0.76,0.045]); sm=plt.cm.ScalarMappable(norm=nz,cmap=cmz)
cb=fig.colorbar(sm,cax=cax_,orientation='horizontal'); cb.set_ticks([-6,-4,-2,0,2]); cb.ax.tick_params(labelsize=6,length=2)
cb.set_label(r'$\log_{10}(\rho/\rho_{\rm unif})$',fontsize=6.5,labelpad=1)
pr=D['azprof']; nn=pr['n']; ks=np.arange(1,nn+1)
from math import comb
tail=[sum(comb(nn,m) for m in range(k,nn+1))/2**nn for k in ks]; pmf=[comb(nn-1,k-1)/2**nn for k in ks]
axs[1].plot(ks,pr['sw'],'o',ms=2.5,color=C1,label='adjacent side (exact)'); axs[1].plot(ks,tail,'-',color=C1,lw=0.8,label=r'$P(\mathrm{Bin}(n,\frac{1}{2})\geq k)$')
axs[1].plot(ks,np.array(pr['ne'])*10,'s',ms=2.5,color=C2,label=r'far side $\times10$'); axs[1].plot(ks,np.array(pmf)*10,'-',color=C2,lw=0.8,label=r'$\binom{n-1}{k-1}2^{-n}\times10$')
axs[1].set_xlabel('position $k$ along the side'); axs[1].set_ylabel('weight'); axs[1].set_xlim(-1,66); axs[1].legend(frameon=False,fontsize=5.5,loc='upper right',handlelength=1.2); axs[1].set_title(f'(b) exact boundary laws, $n={nn}$')
m=np.array(D['azmass']); loglog_series(axs[2],m[:,0],1-m[:,1],C1,'off the adjacent side',-1); loglog_series(axs[2],m[:,0],m[:,3],C2,'interior'); loglog_series(axs[2],m[:,0],m[:,2],C3,'far side')
axs[2].set_xlabel('$n$'); axs[2].set_ylabel('total weight'); axs[2].legend(frameon=False,fontsize=6,loc='lower left'); axs[2].minorticks_off(); axs[2].set_xticks([int(v) for v in m[:,0]][::2]); axs[2].set_xticklabels([str(int(v)) for v in m[:,0]][::2]); axs[2].set_title('(c) boundary-layer trapping')
save(fig,'fig6_aztec')
# ---------- Fig 7 Aztec mu ----------
mu=D['mu']; fig,axs=plt.subplots(1,3,figsize=(7,2.35),gridspec_kw=dict(wspace=0.35))
ax=axs[0]; ax.plot([-1,0,1,0,-1],[0,1,0,-1,0],'k',lw=0.9); th=np.linspace(0,2*np.pi,300); r=1/np.sqrt(2)
ax.fill(r*np.cos(th),r*np.sin(th),color='#dbe9f6'); ax.plot(r*np.cos(th),r*np.sin(th),color=C1,lw=0.8)
k=0.2; P=np.array([-1+k,-k]); Dd=2*(2*k*k-2*k+1); Q=np.array([(2*k*k-1)/Dd,(2*k*k-4*k+1)/Dd])
dv=Q-P; ax.plot([P[0],P[0]+1.7*dv[0]],[P[1],P[1]+1.7*dv[1]],color=C2,lw=1)
ax.plot(*P,'o',color=C2,ms=4); ax.plot(*Q,'o',color=C1,ms=4); ax.text(P[0]-0.05,P[1]-0.16,r'$P(\kappa)$',ha='right'); ax.text(Q[0]+0.06,Q[1]-0.02,r'$Q(\kappa)$')
ax.text(0,0,'liquid',ha='center',va='center'); ax.text(-0.8,-0.02,'frozen',fontsize=6.5,ha='center')
ax.set_aspect('equal'); ax.axis('off'); ax.set_title('(a) tangent continuation')
me=np.array(mu['mu_exact']); mf=np.array(mu['mu_fem'])
axs[1].plot(me[:,0],me[:,1],'-',color='k',lw=0.9,label=r'$\ln\frac{1-\kappa}{\kappa}$ (exact lattice)'); axs[1].plot(mf[:,0],mf[:,1],'o',ms=2.5,color=C1,label='limit shape + tangent')
axs[1].set_xlabel(r'$\kappa$'); axs[1].set_ylabel(r'$\mu(\kappa)\,\delta/2$'); axs[1].legend(frameon=False,fontsize=6); axs[1].set_title('(b) boundary chemical potential')
ce=np.array(mu['cum_exact']); cf=np.array(mu['cum_fem']); cn=np.array(mu['cum_naive'])
axs[2].plot(ce[:,0],ce[:,1],'-',color='k',lw=0.9,label=r'$\ln2-H(\kappa)$'); axs[2].plot(cf[:,0],cf[:,1],'o',ms=2.5,color=C1,label='tangent continuation')
axs[2].plot(cn[:,0],cn[:,1],'s',ms=2.5,color=C2,label='naive penalty FEM')
axs[2].set_xlabel(r'$\kappa$'); axs[2].set_ylabel(r'$\frac{1}{n}[\ln w(\frac{1}{2})-\ln w(\kappa)]$'); axs[2].legend(frameon=False,fontsize=6); axs[2].set_title('(c) integrated rate')
save(fig,'fig7_aztec_mu')
# ---------- Fig 8 house ----------
d=pickle.load(open('fem_house64_0.5.pkl','rb')); Nh=d['N']; g=d['g']; Cc=d['C']; l1=np.abs(g).sum(1)
fig,axs=plt.subplots(1,2,figsize=(7,2.9),gridspec_kw=dict(width_ratios=[1,1.3],wspace=0.3))
ax=axs[0]; col=np.where(l1<0.95,0.0,1.0)
phase=np.full(len(g),0); PH={(1,0):1,(-1,0):2,(0,1):3,(0,-1):4}
for i,(s,t) in enumerate(g):
    if l1[i]>=0.95: phase[i]=PH.get((int(round(s)),int(round(t))),0)
cols=np.array(['#dbe9f6','#f4a582','#92c5de','#d6d6d6','#b8e186'])
ax.scatter(Cc[:,0],Cc[:,1],c=cols[phase],s=1.2,marker='s',linewidths=0)
k=Nh//2; xs=int(0.078*Nh)
ax.plot([0,0,xs],[0,Nh-k,Nh-k+xs],color=C2,lw=1.6); ax.plot(0,0,'o',color=C2,ms=4); ax.plot(xs,Nh-k+xs,'o',color=C2,ms=4)
ax.text(1.5,1.5,'a1',fontsize=7); ax.text(xs+1,Nh-k+xs+1,'optimal defect',fontsize=6.5)
ax.text(Nh*0.5,Nh*0.4,'liquid',ha='center'); ax.text(Nh*0.5,Nh*0.85,'frozen',ha='center')
ax.set_aspect('equal'); ax.axis('off'); ax.set_title(r'(a) house, $k=N/2$: path $\Gamma$')
H2=D['houseR']['k2']['rows']; H4=D['houseR']['k4']['rows']; H2=np.array(H2); H4=np.array(H4)
axs[1].plot(H2[:,0],H2[:,1],'-',color=C2,label=r'prediction $k=N/2$'); axs[1].plot(H4[:,0],H4[:,1],'-',color=C1,label=r'prediction $k=N/4$')
axs[1].axhline(0.213,ls='--',color=C2,lw=0.8,label='exact lattice rate'); axs[1].axhline(0.255,ls='--',color=C1,lw=0.8)
axs[1].axvspan(0.08,0.09,color='gray',alpha=0.2,label='exact optimum')
axs[1].set_xlabel('endpoint on the staircase $x/N$'); axs[1].set_ylabel(r'$R(x)=\frac{1}{N}\int_\Gamma\mu\,ds$'); axs[1].legend(frameon=False,fontsize=6); axs[1].set_title('(b) parameter-free rate prediction')
save(fig,'fig8_house')
# ---------- Fig 9 lozenge ----------
if os.path.exists('lozfig.json'):
  L=json.load(open('lozfig.json'))
  fig,axs=plt.subplots(1,2,figsize=(7,2.6),gridspec_kw=dict(wspace=0.3))
  for (key,col,lab) in (('1,1,1',C1,'$1{:}1{:}1$'),('1,1.5,0.75',C2,'$1{:}1.5{:}0.75$')):
      rows=L[key]['rows']; kk=[r['kap'] for r in rows]
      axs[0].plot(kk,[r['rule'] for r in rows],'-',color=col,label='rule, hexagon '+lab); axs[0].plot(kk,[r['tau1'] for r in rows],':',color=col,label=r'$\tau=1$ only')
      axs[0].plot(kk,[r['extrap'] for r in rows],'o',color=col,ms=3.5,label='exact (extrap.)')
  axs[0].set_xlabel(r'dent positions $\kappa=\lambda$'); axs[0].set_ylabel('rate $\\frac{1}{n}\\ln w$'); axs[0].legend(frameon=False,fontsize=5.8); axs[0].set_title('(a) dent pair at a corner')
  hexs=['1:1:1','1:2:1.5','2:1:1.5','1.5:1:2','0.6:1.4:1']; need=[0.52325,0.70115,0.70115,0.47414,0.44271]
  da=[0.52289,1.18115,0.47394,0.70082,0.85032]; db=[0.52289,0.47394,1.18115,1.18115,0.27931]; dc=[0.52289,0.70082,0.70082,0.47394,0.44240]
  x=np.arange(5); axs[1].bar(x-0.27,da,0.18,color='#bbbbbb',label=r'$\partial_a\ln Z/n$'); axs[1].bar(x-0.09,db,0.18,color='#888888',label=r'$\partial_b\ln Z/n$')
  axs[1].bar(x+0.09,dc,0.18,color=C1,label=r'$\partial_c\ln Z/n$'); axs[1].plot(x+0.27,need,'*',color=C2,ms=7,label='required layer cost')
  axs[1].set_xticks(x); axs[1].set_xticklabels(hexs,fontsize=6,rotation=20); axs[1].set_xlabel('hexagon $a{:}b{:}c$'); axs[1].legend(frameon=False,fontsize=6); axs[1].set_title('(b) layer cost $=\\partial_c\\ln Z$')
  save(fig,'fig9_lozenge')
# ---------- Fig 10 triangular ----------
fig,axs=plt.subplots(1,2,figsize=(7,2.5),gridspec_kw=dict(width_ratios=[1,1.2],wspace=0.3))
T=np.array(D['tri13']); ax=axs[0]; vmax=2.0
for i in range(T.shape[0]):
    for j in range(T.shape[1]):
        x0=i+0.5*j; y0=j*np.sqrt(3)/2
        ax.add_patch(Polygon([[x0,y0],[x0+1,y0],[x0+1.5,y0+np.sqrt(3)/2],[x0+0.5,y0+np.sqrt(3)/2]],fc=cmap(min(T[i,j]/vmax,1)),ec='white',lw=0.2))
ax.set_xlim(0,20); ax.set_ylim(0,12); ax.set_aspect('equal'); ax.axis('off'); ax.set_title('(a) triangular lattice, $13\\times13$ rhombus')
for key,col,lab in (('t60',C2,r'triangular, $60^\circ$ corner'),('t120',C1,r'triangular, $120^\circ$ corner'),('sq',C3,'square, odd board')):
    a=np.array(D['triratio'][key]); axs[1].loglog(a[:,0],a[:,1],'o-',color=col,ms=3,label=lab)
axs[1].set_xlabel('size $N$'); axs[1].set_ylabel('corner/centre ratio'); axs[1].legend(frameon=False,fontsize=6.5); axs[1].set_title('(b) no long-range corner effect without charge')
save(fig,'fig10_triangular')
