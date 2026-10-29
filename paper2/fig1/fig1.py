import json, sys, numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle, FancyArrowPatch
from matplotlib.collections import LineCollection
REPO='/tmp/claude-0/-home-claude/76b76c06-8fe9-5f68-b413-27ba31b71e09/scratchpad/repo/domino-defects-code/paper2'
sys.path.insert(0,REPO)
from honey_defect import flakeA, D
plt.rcParams.update({'font.size':10,'font.family':'DejaVu Sans','mathtext.fontset':'dejavusans'})
INK='#1d2733'; MUTE='#6b7680'; BLUE='#2b6cb0'; ORANGE='#dd6b20'; RED='#c0392b'; H1='#dbe7f3'; H2='#f6e6d6'
fig=plt.figure(figsize=(7.0,4.8))
axa=fig.add_axes([0.02,0.40,0.40,0.535]); axb=fig.add_axes([0.45,0.40,0.53,0.535])
axc=[fig.add_axes([0.02+0.33*k,0.165,0.10,0.135]) for k in range(3)]
for ax in [axa,axb]+axc: ax.set_axis_off()
# ---------- (a) the problem
s=json.load(open('sample.json')); L=8
for a,b in s['dominoes']:
    (r1,c1),(r2,c2)=a,b; x0,y0=min(c1,c2),min(r1,r2); w=1+abs(c1-c2); h=1+abs(r1-r2)
    axa.add_patch(FancyBboxPatch((x0+0.09,y0+0.09),w-0.18,h-0.18,boxstyle='round,pad=0,rounding_size=0.12',
                  fc=H1 if w==2 else H2,ec=INK,lw=0.8))
for r in range(L):
    for c in range(L):
        if (r,c)!=(0,0) and (r+c)%2==0: axa.add_patch(Rectangle((c+0.42,r+0.42),0.16,0.16,fc='#9aa5b1',ec='none',zorder=3))
mr,mc=s['mono']; axa.add_patch(Rectangle((mc+0.06,mr+0.06),0.88,0.88,fc='none',ec=RED,lw=1.8,zorder=4))
axa.text(mc+0.5,mr+0.5,'?',ha='center',va='center',color=RED,fontsize=13,fontweight='bold',zorder=5)
axa.add_patch(Rectangle((0,0),1,1,fc='white',ec=INK,lw=0.8,hatch='////',zorder=2))
axa.plot([0,L,L,0,0],[0,0,L,L,0],color=INK,lw=1.2)
axa.set_xlim(-0.3,L+0.3); axa.set_ylim(-1.75,L+0.3); axa.set_aspect('equal')
axa.text(L/2,-0.45,'Hatched corner removed: one more white than\nblack cell, so one cell must stay uncovered.\nWhere does it sit?',ha='center',va='top',fontsize=7.5,color=INK)
# ---------- (b) the answer: exact density + electrostatics
d=json.load(open('density.json')); Lb=d['L']; W=np.array(d['W']); p=np.array(d['p'])
img=np.full((Lb,Lb),np.nan); 
for (r,c),v in zip(W,p): img[r,c]=np.log10(v)
for (r,c) in W:                      # fill black cells with neighbour mean for a smooth picture
    pass
fill=img.copy()
for r in range(Lb):
    for c in range(Lb):
        if np.isnan(img[r,c]) and (r,c)!=(0,0):
            nbv=[img[x,y] for x,y in [(r+1,c),(r-1,c),(r,c+1),(r,c-1)] if 0<=x<Lb and 0<=y<Lb and not np.isnan(img[x,y])]
            fill[r,c]=np.mean(nbv)
im=axb.imshow(fill,origin='lower',extent=(0,Lb,0,Lb),cmap='magma_r',interpolation='nearest')
axb.add_patch(Rectangle((0,0),1,1,fc='white',ec=INK,lw=0.6,hatch='////'))
axb.plot([0,Lb,Lb,0,0],[0,0,Lb,Lb,0],color=INK,lw=1.2)
def charge(x,y,q,col):
    axb.add_patch(Circle((x,y),1.55,fc='white',ec=col,lw=1.3,zorder=5)); axb.text(x,y,q,ha='center',va='center',fontsize=7.6,color=col,zorder=6,fontweight='bold')
charge(-1.2,-1.2,'+¾',RED); charge(Lb+1.2,-1.2,'+¼',RED); charge(-1.2,Lb+1.2,'+¼',RED); charge(Lb+1.2,Lb+1.2,'−¼',BLUE)
# image across the right edge
ex,ey=Lb-5.5,Lb*0.55
axb.add_patch(Circle((ex,ey),0.9,fc=RED,ec='white',lw=0.8,zorder=6)); axb.text(ex,ey,'e',color='white',ha='center',va='center',fontsize=7,zorder=7)
axb.add_patch(Circle((2*Lb-ex,ey),0.9,fc='white',ec=RED,lw=1.0,ls='--',zorder=6)); axb.text(2*Lb-ex,ey,'e',color=RED,ha='center',va='center',fontsize=7,zorder=7)
axb.annotate('',xy=(2*Lb-ex-1.0,ey),xytext=(ex+1.0,ey),arrowprops=dict(arrowstyle='<->',color=RED,lw=0.8,ls=':'),zorder=6)
axb.text(Lb+1.5,ey+2.2,'edge = mirror,\nimage of\nsame sign',fontsize=7,color=MUTE,ha='left',va='bottom')
axb.set_xlim(-3.2,Lb+9.5); axb.set_ylim(-5.0,Lb+3.2); axb.set_aspect('equal')
cax=fig.add_axes([0.60,0.41,0.16,0.016]); cb=fig.colorbar(im,cax=cax,orientation='horizontal'); cb.set_ticks([]); cb.outline.set_linewidth(0.4)
fig.text(0.595,0.418,'low',fontsize=7,color=MUTE,ha='right',va='center'); fig.text(0.765,0.418,'high',fontsize=7,color=MUTE,ha='left',va='center')
fig.text(0.68,0.39,r'exact probability of the uncovered cell',fontsize=7,color=MUTE,ha='center')
# ---------- (c) three knobs
k0,k1,k2=axc
# stiffness: plaquette with two parallel dimers
for (x,y) in [(0,0),(0,1)]:
    k0.add_patch(FancyBboxPatch((x+0.08,y+0.08),1.84,0.84,boxstyle='round,pad=0,rounding_size=0.12',fc=H1,ec=INK,lw=0.8))
k0.add_patch(Rectangle((0.5,0.5),1,1,fc='none',ec=ORANGE,lw=1.2,ls='--'))
k0.text(1.0,1.0,'v',color=ORANGE,ha='center',va='center',fontsize=9,fontweight='bold')
k0.set_xlim(-0.2,2.2); k0.set_ylim(-0.2,2.2); k0.set_aspect('equal')
# bipartiteness: square lattice with diagonals
xs,ys=np.meshgrid(range(3),range(3)); segs=[]
for x in range(3):
    for y in range(3):
        if x<2: segs.append([(x,y),(x+1,y)])
        if y<2: segs.append([(x,y),(x,y+1)])
k1.add_collection(LineCollection(segs,colors=INK,lw=0.9))
k1.add_collection(LineCollection([[(x,y),(x+1,y+1)] for x in range(2) for y in range(2)],colors=ORANGE,lw=1.1,linestyles='--'))
k1.scatter(xs,ys,s=10,c=INK,zorder=3); k1.text(1.55,0.25,'t',color=ORANGE,fontsize=9,fontweight='bold')
k1.set_xlim(-0.3,2.3); k1.set_ylim(-0.3,2.3); k1.set_aspect('equal')
# lattice: honeycomb armchair flake corner region, boundary coloured by edge class
A,B=flakeA(4.0); key=lambda q:(round(q[0],4),round(q[1],4)); Bs={key(b) for b in B}
k2.add_collection(LineCollection([[a,a+d] for a in A for d in D if key(a+d) in Bs],colors='#9aa5b1',lw=0.7))
for P in (A,B):
    pr=np.array([[x@np.array([np.cos(np.radians(t)),np.sin(np.radians(t))]) for t in range(0,360,60)] for x in P])
    mx=pr.max(1); th=np.arange(0,360,60)[pr.argmax(1)]; bd=mx>mx.max()-1.2
    k2.scatter(P[bd,0],P[bd,1],s=7,c=np.where(th[bd]%120==0,BLUE,ORANGE),zorder=3)
k2.set_aspect('equal'); k2.autoscale_view()
txt=[('Stiffness $g$  (dimer interactions)',r'charges fixed, coupling $\propto g$:'+'\n'+r'every exponent rescales by $2g$   (Sec. II)'),
     ('Bipartiteness  (diagonal bonds $t$)','no height field: the force is'+'\n'+r'screened beyond $\xi\sim1/t$   (Sec. III)'),
     ('Lattice  (honeycomb)','corner charge read from the boundary'+'\n'+r'height, shifted by $\pm\frac{1}{4}$ by edge class   (Sec. IV)')]
for k,(h,b) in enumerate(txt):
    x=0.135+0.33*k
    fig.text(x,0.285,h.replace('  (','\n('),fontsize=8.2,fontweight='bold',color=INK,va='top',linespacing=1.2)
    fig.text(0.02+0.33*k,0.145,b.replace('   (Sec.','\n(Sec.'),fontsize=7.4,color=INK,va='top',linespacing=1.3)
for x,y,lab,t in [(0.02,0.985,'(a)','The question'),(0.45,0.985,'(b)','The answer: two-dimensional electrostatics'),(0.02,0.335,'(c)','What controls it')]:
    fig.text(x,y,lab,fontsize=9.5,fontweight='bold',va='top',color=INK); fig.text(x+0.04,y,t,fontsize=9.5,va='top',color=INK)
fig.add_artist(plt.Line2D([0.02,0.98],[0.355,0.355],color='#d5dbe1',lw=0.8))
fig.savefig('fig00_overview.pdf'); fig.savefig('fig00_overview.png',dpi=200)
