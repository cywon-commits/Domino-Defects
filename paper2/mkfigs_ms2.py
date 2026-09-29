import json, numpy as np, sys, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
sys.path.insert(0,'.')
from honeyflake import flake, probs, D
from honey_defect import flakeA, Flake, bonds
from honey_triangle import triA
from honey_para import paraA
from honey_corner import boundary
plt.rcParams.update({'font.size':13,'axes.spines.top':False,'axes.spines.right':False})
F='manuscript/figs/'; C1,C2,C3='#1f3a5f','#e07b39','#8a3a12'
def draw(ax,A,B,c0,title,normals):
    key=lambda p:(round(p[0],4),round(p[1],4)); Bs={key(b) for b in B}
    segs=[[a,a+d] for a in A for d in D if key(a+d) in Bs]
    ax.add_collection(LineCollection(segs,colors='#9aa5b1',lw=0.8))
    for P,col in [(A,C1),(B,'white')]:
        ax.scatter(P[:,0],P[:,1],s=10,c=col,edgecolors=C1,lw=0.5,zorder=3)
    # boundary vertices coloured by edge class
    for P in (A,B):
        X=P-c0; pr=np.array([[x@np.array([np.cos(np.radians(t)),np.sin(np.radians(t))]) for t in normals] for x in X])
        mx=pr.max(1); th=np.array(normals)[pr.argmax(1)]; bd=mx>mx.max()-1.2
        cls=np.where(th%120==0,'#2b6cb0','#dd6b20')
        ax.scatter(P[bd,0],P[bd,1],s=26,c=cls[bd],zorder=4)
    ax.set_aspect('equal'); ax.axis('off'); ax.set_title(title,fontsize=12)
fig,ax=plt.subplots(1,3,figsize=(12,4.3))
A,B=flakeA(7); draw(ax[0],A,B,np.zeros(2),'hexagon: 120° corners',[0,60,120,180,240,300])
A,B=triA(9); draw(ax[1],A,B,np.zeros(2),'triangle: 60° corners',[0,120,240])
A,B=paraA(6,np.zeros(2)); draw(ax[2],A,B,np.zeros(2),'parallelogram: 60° and 120°',[0,60,180,240])
ax[2].scatter([],[],c='#2b6cb0',label='edge class + (normal 0°,120°,240°)'); ax[2].scatter([],[],c='#dd6b20',label='edge class − (60°,180°,300°)')
ax[2].legend(loc='lower center',bbox_to_anchor=(-0.7,-0.18),ncol=2,frameon=False,fontsize=11)
plt.tight_layout(); plt.savefig(F+'fig01b_honey.pdf'); plt.close()
# flatness maps
fig,ax=plt.subplots(1,2,figsize=(10,4.8))
for k,(kind,R) in enumerate([('armchair',20),('zigzag',20)]):
    A,B=flake(R,kind); E,P=probs(A,B); pos=np.array([(A[e[0]]+B[e[1]])/2 for e in E])
    sc=ax[k].scatter(pos[:,0],pos[:,1],c=np.abs(P-1/3),s=9,cmap='magma_r',vmin=0,vmax=0.5)
    ax[k].set_aspect('equal'); ax[k].axis('off'); ax[k].set_title('%s flake, R = %d'%(kind,R))
fig.colorbar(sc,ax=ax,shrink=0.8,label=r'$|P_{\rm bond}-1/3|$'); plt.savefig(F+'fig06_flatness.pdf',bbox_inches='tight'); plt.close()
# boundary heights
fig,ax=plt.subplots(figsize=(8.5,4.2))
for (lab,R,col) in [('type I, R = 28',28,C1),('type II, R = 27',27,C2)]:
    A,B=flakeA(R); X,H,tot=boundary(A,B); ang=np.degrees(np.arctan2(X[:,1],X[:,0]))%360; o=np.argsort(ang)
    ax.plot(ang[o],H[o],'.',ms=3,color=col,label=lab)
for c in range(30,360,60): ax.axvline(c,color='#bbbbbb',lw=0.8,ls=':')
ax.set_xlabel('polar angle of boundary step (deg)'); ax.set_ylabel('boundary height $h$ (midpoints)'); ax.legend(frameon=False)
plt.tight_layout(); plt.savefig(F+'fig07_boundary_height.pdf'); plt.close()
# density maps
fig,ax=plt.subplots(1,2,figsize=(10,4.8))
for k,(lab,R) in enumerate([('type I (Temperley), R = 20',20),('type II, R = 22',22)]):
    Fl=Flake(R); lz=np.array([Fl.logZ_minus(v) for v in range(len(Fl.A))]); p=np.exp(lz-lz.max()); p/=p.sum()
    sc=ax[k].scatter(Fl.A[:,0],Fl.A[:,1],c=np.log10(p),s=16,cmap='viridis'); ax[k].set_aspect('equal'); ax[k].axis('off'); ax[k].set_title(lab)
    fig.colorbar(sc,ax=ax[k],shrink=0.75,label=r'$\log_{10}\rho$')
plt.savefig(F+'fig08_density.pdf',bbox_inches='tight'); plt.close()
# exponents
D_='../data/'; hs=json.load(open(D_+'paper2_honey_scaling.json')); tr=np.array(json.load(open(D_+'paper2_honey_triangle.json'))); pa=json.load(open(D_+'paper2_honey_para.json'))
fig,ax=plt.subplots(1,3,figsize=(14,4.3))
def panel(a,R,y,lab,g,col,mk):
    a.loglog(R,y,mk,color=col,ms=6,label=lab+r'  (pred. $\gamma$=%.2f)'%g,mfc='none',mew=1.6)
    RR=np.array([R.min(),R.max()]); a.loglog(RR,y[-1]*(RR/R[-1])**(-g),'-',color=col,lw=1)
for typ,col in [('I',C1),('II',C2)]:
    a=np.array(hs[typ]); gp={'I':(0.25,0.25),'II':(0.75,-0.25)}[typ]
    panel(ax[0],a[:,0],a[:,1],'type %s, 30°'%typ,gp[0],col,'o'); panel(ax[0],a[:,0],a[:,2],'type %s, 90°'%typ,gp[1],col,'s')
panel(ax[1],tr[:,0],tr[:,2],'60° corner',0.25,C1,'o')
for fam,gs,col in [('F1',(0.75,0.25,-0.25,0.25),C1),('F2',(0.25,0.25,0.25,0.25),C2)]:
    a=np.array(pa[fam]['data'])
    for j,(c,mk) in enumerate([(30,'o'),(120,'s'),(210,'^'),(300,'D')]):
        if fam=='F2' and j>0: continue
        panel(ax[2],a[:,0],a[:,j+1],'%s %d°'%(fam,c),gs[j],col if j<2 else C3,mk)
for a,t in zip(ax,['hexagons (120° corners)','triangles (60° corners)','parallelograms']):
    a.set_xlabel('R'); a.set_ylabel('corner / centre'); a.set_title(t,fontsize=12); a.legend(frameon=False,fontsize=8.5)
plt.tight_layout(); plt.savefig(F+'fig09_exponents.pdf'); plt.close()
# predicted vs measured
def gfit(R,y,n=5): return -np.polyfit(np.log(R[-n:]),np.log(y[-n:]),1)[0]
pts=[]
a=np.array(hs['I']); pts+= [(0.25,0.5,gfit(a[:,0],a[:,1]),'hex I'),(0.25,0.0,gfit(a[:,0],a[:,2]),'hex I')]
a=np.array(hs['II']); pts+= [(0.75,1.0,gfit(a[:,0],a[:,1]),'hex II'),(-0.25,-0.5,gfit(a[:,0],a[:,2]),'hex II')]
pts+= [(0.25,0.25,gfit(tr[:,0],tr[:,2],6),'triangle')]
for fam,gs,gn in [('F1',(0.75,0.25,-0.25,0.25),(1,0.25,-0.5,0.25)),('F2',(0.25,)*4,(0.5,0.25,0,0.25))]:
    a=np.array(pa[fam]['data'])
    for j in range(4): pts.append((gs[j],gn[j],gfit(a[:,0],a[:,j+1]),'para '+fam))
pts=np.array([(p[0],p[1],p[2]) for p in pts])
fig,ax=plt.subplots(figsize=(5.6,5.2))
ax.plot([-0.7,1.2],[-0.7,1.2],'k:',lw=1)
ax.plot(pts[:,1],pts[:,2],'x',color='#999999',ms=8,mew=1.6,label='naive rule $q_c=-J/3$')
ax.plot(pts[:,0],pts[:,2],'o',color=C1,ms=7,label='edge-class rule $q_c=-(J+\\sigma/2)/3$')
ax.set_xlabel(r'predicted $\gamma$'); ax.set_ylabel(r'measured $\gamma$ (largest 5 R)'); ax.legend(frameon=False,fontsize=10,loc='upper left')
plt.tight_layout(); plt.savefig(F+'fig10_pred_vs_meas.pdf'); plt.close()
json.dump(pts.tolist(),open(D_+'paper2_pred_vs_meas.json','w'))
print('pred/naive/meas:\n',np.round(pts,3))
# J_eff and strip
fams={'J=0 (pred +1/2)':([12,17,22,25,27,30],[0.32,0.69,0.37,0.42,0.42,0.46],0.5),
      'J=-1 (pred -1/2)':([11,16,21,24,26,29],[-0.24,-0.35,-0.61,-0.61,-0.42,-0.46],-0.5),
      'J=-2 (pred -3/2)':([15,18,20,23,28],[-1.32,-1.33,-1.09,-1.39,-1.24],-1.5)}
fig,ax=plt.subplots(1,2,figsize=(11,4.2))
for (lab,(R,v,p)),col in zip(fams.items(),[C1,C2,C3]):
    ax[0].plot(1/np.array(R),v,'o',color=col,label=lab); ax[0].axhline(p,color=col,ls='--',lw=1); ax[0].axhline(round(p-0.5 if p>0 else p+0.5),color=col,ls=':',lw=1)
ax[0].set_xlabel('1/R'); ax[0].set_ylabel(r'$J_{\rm eff}$ at the 30° corner'); ax[0].legend(frameon=False,fontsize=9); ax[0].set_title('dashed: edge-class rule, dotted: naive $J$',fontsize=10)
st=json.load(open(D_+'paper2_honey_strip2.json'))
for W,L,dc,d in st: ax[1].plot(W+(0.3 if dc==5 else 0),d,'o' if dc==3 else 's',color=C1 if dc==3 else C2,ms=7)
ax[1].axhline(-0.5,color='k',ls='--',lw=1); ax[1].set_xlabel('strip half-width W'); ax[1].set_ylabel(r'$\delta_+-\delta_-$'); ax[1].set_ylim(-0.8,-0.3)
ax[1].plot([],[],'o',color=C1,label='edge layer cut 3'); ax[1].plot([],[],'s',color=C2,label='edge layer cut 5'); ax[1].legend(frameon=False,fontsize=9)
plt.tight_layout(); plt.savefig(F+'fig11_Jeff.pdf'); plt.close()
# image test
im=json.load(open(D_+'paper2_honey_image.json'))
fig,ax=plt.subplots(1,2,figsize=(11,4.2))
for typ,col in [('I',C1),('II',C2)]:
    a=np.array(im[typ]); ax[0].loglog(a[:,0],a[:,1],'o',color=col,label='type '+typ,mfc='none',mew=1.5)
RR=np.array([9,38.]); ax[0].loglog(RR,0.83*(RR/9)**-0.25,'k--',lw=1,label=r'$R^{-1/4}$ (same sign)'); ax[0].loglog(RR,0.83*(RR/9)**0.25,'k:',lw=1,label=r'$R^{+1/4}$ (opposite)')
ax[0].set_xlabel('R'); ax[0].set_ylabel('edge midpoint / centre'); ax[0].legend(frameon=False,fontsize=9)
for typ,col in [('I',C1),('II',C2)]:
    key=[k for k in im if k.startswith('profile_'+typ+'_')][0]; d,r=im[key]; ax[1].loglog(d,r,'.',color=col,label='type '+typ)
dd=np.array([2,10.]); ax[1].loglog(dd,0.6*(dd/2)**0.25,'k--',lw=1,label=r'$d^{+1/4}$')
ax[1].set_xlabel('distance d from the edge midpoint'); ax[1].set_ylabel(r'$\rho(d)/\rho({\rm centre})$'); ax[1].legend(frameon=False,fontsize=9)
plt.tight_layout(); plt.savefig(F+'fig12_image.pdf'); plt.close()
