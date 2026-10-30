import json, numpy as np, sys, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0,'.'); from sc import density, inv
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
d=json.load(open('exact_28_27_36_38.json')); C1,C2,GR='#1f3a5f','#dd6b20','#8a8f96'
fig,ax=plt.subplots(1,3,figsize=(13,3.9))
for k,(R,typ,qe,qn) in enumerate([(36,'I',[-1/6]*6,[0,-1/3]*3),(38,'II',[1/6,-1/2]*3,[1/3,-2/3]*3)]):
    A=np.array(d[str(R)]['A']); lz=np.array(d[str(R)]['lz']); Z=A[:,0]+1j*A[:,1]
    hexr=np.max(np.stack([A@np.array([np.cos(np.radians(60*j)),np.sin(np.radians(60*j))]) for j in range(6)]),axis=0)
    pre=inv(Z,R+0.5); le,_=density(Z,R+0.5,qe,pre=pre); ln,_=density(Z,R+0.5,qn,pre=pre)
    m=hexr<R-3; c0=np.mean(lz[m]-le[m]); c1=np.mean(lz[m]-ln[m])
    for lab,ang,mk in [('to 30° corner',30,'o'),('to 90° corner',90,'s'),('to edge midpoint',0,'^')]:
        u=np.array([np.cos(np.radians(ang)),np.sin(np.radians(ang))]); s=A@u; perp=np.abs(A@np.array([-u[1],u[0]]))
        sel=np.where(perp<0.9)[0]; sel=sel[np.argsort(s[sel])]; sel=sel[s[sel]>0]
        x=s[sel]/(R/np.cos(np.pi/6) if ang!=0 else R)
        ax[k].plot(x,np.exp(lz[sel]-c0),mk,ms=4,color=C1 if ang==30 else (C2 if ang==90 else GR),mfc='none',mew=1.2,label='exact, '+lab)
        ax[k].plot(x,np.exp(le[sel]),'-',color=C1 if ang==30 else (C2 if ang==90 else GR),lw=1.4)
        ax[k].plot(x,np.exp(ln[sel]+c1-c0),':',color=C1 if ang==30 else (C2 if ang==90 else GR),lw=1.2)
    ax[k].set_xlabel('distance from centre / distance to boundary'); ax[k].set_ylabel(r'$\rho$ (common normalisation)')
    ax[k].set_title('type %s hexagon, R = %d'%(typ,R),fontsize=11)
ax[0].legend(frameon=False,fontsize=8,loc='lower left')
ax[1].plot([],[],'k-',lw=1.4,label='edge-class rule'); ax[1].plot([],[],'k:',lw=1.2,label='naive rule'); ax[1].legend(frameon=False,fontsize=9,loc='upper left')
c=json.load(open('compare.json')); xs=[];
for j,(R,typ) in enumerate([(28,'I'),(36,'I'),(27,'II'),(38,'II')]):
    rows=[r for r in c[str(R)] if r[0]==0.5 and r[2]==6]
    for r in rows: ax[2].bar(j+(-0.18 if r[1]=='edge-class' else 0.18),r[3],width=0.34,color=C1 if r[1]=='edge-class' else '#c9cfd6')
ax[2].set_xticks(range(4)); ax[2].set_xticklabels(['I, 28','I, 36','II, 27','II, 38']); ax[2].set_ylabel(r'rms deviation of $\ln\rho$ (interior)')
from matplotlib.patches import Patch
ax[2].legend(handles=[Patch(color=C1,label='edge-class rule'),Patch(color='#c9cfd6',label='naive rule')],frameon=False,fontsize=9,loc='upper left')
ax[2].set_ylim(0,0.11)
ax[2].set_title('no adjustable parameter',fontsize=11)
plt.tight_layout(); plt.savefig('fig_hexdens.pdf'); plt.savefig('fig_hexdens.png',dpi=110)
