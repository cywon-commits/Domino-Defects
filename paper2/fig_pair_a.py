import json,numpy as np,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':8.5,'axes.spines.top':False,'axes.spines.right':False})
cols=['#1f3a5f','#3a6ea5','#7fa7d1','#dd6b20','#b5451b','#8a3a12']
p=json.load(open('../data/paper2_pair_bulk.json'))
fig,ax=plt.subplots(figsize=(2.9,2.6))
for t,cl in zip(['0','0.1','0.2','0.5','1.0'],[cols[0],cols[2],cols[3],cols[4],cols[5]]):
    a=np.array(p[t]); m=a[:,1]>0; ax.plot(a[m,0],a[m,1],'o-',ms=3,lw=1.2,color=cl,label='t = %s'%t)
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('separation r'); ax.set_ylabel('C(r)'); ax.set_title('(a) pair function, 64×64',fontsize=8)
ax.legend(frameon=False,fontsize=7); plt.tight_layout(); plt.savefig('manuscript/figs/fig03a_pair.pdf')
