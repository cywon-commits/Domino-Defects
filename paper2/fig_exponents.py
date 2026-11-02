import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, FixedLocator, FixedFormatter
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
C1,C2,C3='#1f3a5f','#dd6b20','#2a9d8f'; D_='../data/'
hs=json.load(open(D_+'paper2_honey_scaling.json')); tr=np.array(json.load(open(D_+'paper2_honey_triangle.json'))); pa=json.load(open(D_+'paper2_honey_para.json'))
fig,ax=plt.subplots(1,3,figsize=(13,4.0))
def panel(a,R,y,lab,g,col,mk):
    a.loglog(R,y,mk,color=col,ms=5.5,label=lab+r' (pred. %.2f)'%g,mfc='none',mew=1.4)
    RR=np.array([R.min(),R.max()]); a.loglog(RR,y[-1]*(RR/R[-1])**(-g),'-',color=col,lw=1)
for typ,col in [('I',C1),('II',C2)]:
    a=np.array(hs[typ]); gp={'I':(0.25,0.25),'II':(0.75,-0.25)}[typ]
    panel(ax[0],a[:,0],a[:,1],'type %s, 30°'%typ,gp[0],col,'o')
    if typ=='II': panel(ax[0],a[:,0],a[:,2],'type %s, 90°'%typ,gp[1],col,'s')
panel(ax[1],tr[:,0],tr[:,2],'60° corner',0.25,C1,'o')
for fam,gs,col in [('F1',(0.75,0.25,-0.25,0.25),C1),('F2',(0.25,0.25,0.25,0.25),C2)]:
    a=np.array(pa[fam]['data'])
    for j,(c,mk) in enumerate([(30,'o'),(120,'s'),(210,'^'),(300,'D')]):
        if fam=='F2' and j>0: continue
        if fam=='F1' and j==3: continue
        panel(ax[2],a[:,0],a[:,j+1],'%s, %d° corner'%(fam,c),gs[j],col if j<2 else C3,mk)
titles=['(a) hexagons: 120° corners','(b) triangles: 60° corners','(c) parallelograms']
for a,t in zip(ax,titles):
    a.set_xlabel('R'); a.set_ylabel('corner / centre'); a.set_title(t,fontsize=11)
    a.xaxis.set_minor_formatter(NullFormatter()); a.xaxis.set_major_locator(FixedLocator([10,20,40])); a.xaxis.set_major_formatter(FixedFormatter(['10','20','40']))
    a.yaxis.set_minor_formatter(NullFormatter())
ax[0].legend(frameon=False,fontsize=8,loc='lower left'); ax[1].legend(frameon=False,fontsize=8,loc='upper right')
ax[2].set_ylim(0.12,14); ax[2].legend(frameon=False,fontsize=8,loc='upper left',ncol=2)
from matplotlib.ticker import LogFormatterSciNotation
ax[1].yaxis.set_minor_formatter(matplotlib.ticker.FormatStrFormatter('%.1f'))
plt.tight_layout(); plt.savefig('manuscript/figs/fig09_exponents.pdf'); plt.savefig('fig_exponents.png',dpi=110)
