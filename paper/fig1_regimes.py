# Fig. 1: two regimes. (a) straight edges: power-law density (64x64, corner a1 removed); (b) staircase edges: Aztec diamond, tip removed.
import json, os, sys, numpy as np
sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','lattice'))
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon, Circle
from matplotlib.colors import Normalize
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','data'))
plt.rcParams.update({'font.family':'serif','font.size':8,'axes.titlesize':8,'mathtext.fontset':'cm'})
C2='#c0562b'; D=json.load(open('figdata.json'))
cmz=plt.get_cmap('YlGnBu'); nz=Normalize(-6,2)
fig,axs=plt.subplots(1,2,figsize=(6.4,3.2))
M=np.log10(np.maximum(np.array(D['single64']['data'],float),1e-12)); N=M.shape[0]
ax=axs[0]; ax.imshow(M,origin='lower',cmap=cmz,norm=nz,interpolation='nearest',extent=(0,N,0,N))
ax.add_patch(Rectangle((-1.5,-1.5),2,2,fc=C2,ec='k',lw=0.4,clip_on=False))
ax.text(1.5,-4,'removed corner',color=C2,fontsize=6.5,va='center')
ax.set_xlim(-3,N+1); ax.set_ylim(-6,N+1); ax.set_aspect('equal'); ax.axis('off')
ax.set_title('(a) straight edges ($\\mu=0$): power laws,\n corner charges and images')
az=D['az_single']; n=az['n']; ax=axs[1]
from kreg import aztec as _az
for r,cc,lv in az['log10']:
    x,y=cc+0.5,r+0.5; ax.add_patch(Polygon([(x-1,y),(x,y+1),(x+1,y),(x,y-1)],closed=True,fc=cmz(nz(lv)),ec='none'))
_reg=set(tuple(map(int,v)) for v in _az(n))
for r in range(-2,2*n+2):
    for cc in range(-2,2*n+2):
        if (r,cc) not in _reg: ax.add_patch(Rectangle((cc,r),1,1,fc='white',ec='none',zorder=3))
th=np.linspace(0,2*np.pi,200); ax.plot(n+n/np.sqrt(2)*np.cos(th),n+n/np.sqrt(2)*np.sin(th),'--',c=C2,lw=0.7,zorder=4)
ax.add_patch(Rectangle((0,n),1,1,fc=C2,ec='k',lw=0.4,zorder=6))
ax.text(-1,n-4,'removed\ntip',color=C2,fontsize=6.5,ha='center',va='top',zorder=7)
ax.set_xlim(-7,2*n+1); ax.set_ylim(-1,2*n+1); ax.set_aspect('equal'); ax.axis('off')
ax.set_title('(b) staircase edges ($\\mu\\neq0$): exponential\n trapping in a boundary layer')
cax=fig.add_axes([0.3,0.06,0.4,0.025]); sm=plt.cm.ScalarMappable(norm=nz,cmap=cmz)
cb=fig.colorbar(sm,cax=cax,orientation='horizontal'); cb.set_ticks([-6,-4,-2,0,2]); cb.set_label(r'$\log_{10}(\rho/\rho_{\rm unif})$',fontsize=7)
fig.subplots_adjust(left=0.02,right=0.98,top=0.86,bottom=0.12,wspace=0.08)
fig.savefig('../paper/figs/fig1_regimes.pdf'); fig.savefig('/tmp/fig1_regimes.png',dpi=130)
