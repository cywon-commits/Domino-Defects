import numpy as np, json, os, glob
def load(fn,L): h=np.fromfile(fn).reshape(10,L*L); return h
def corner_ratio(h,L,r=2):
    x,y=np.meshgrid(range(L),range(L),indexing='xy'); x=x.ravel(); y=y.ravel(); b=((x+y)%2==0)
    cm=np.zeros(L*L,bool)
    for cx,cy in [(0,0),(0,L-1),(L-1,0),(L-1,L-1)]: cm|=b&(abs(x-cx)+abs(y-cy)<=r)
    ce=b&(abs(x-L//2)+abs(y-L//2)<=r)
    return (h[:,cm].sum(1)/cm.sum())/(h[:,ce].sum(1)/ce.sum())
def near(h,L,site,R=3.01):
    x=np.arange(L*L)%L; y=np.arange(L*L)//L; sx,sy=site%L,site//L
    m=np.hypot(x-sx,y-sy)<=R; return h[:,m].sum(1)/h.sum(1)
def jk(vals_by_L,Ls):  # fit slope with jackknife over blocks
    Y=np.log(np.array(vals_by_L))  # shape nL x 10
    full=np.polyfit(np.log(Ls),np.log(np.array(vals_by_L).mean(1)),1)[0]
    s=[np.polyfit(np.log(Ls),np.log(np.delete(np.array(vals_by_L),i,1).mean(1)),1)[0] for i in range(10)]
    return full,np.sqrt(9/10*np.sum((np.array(s)-np.mean(s))**2))
out={'corner':{},'vac':{}}
summ=json.load(open('../../data/paper2_worm_summary.json')); gc=dict(zip([float(j) for j in summ['J']],summ['bulk']))
for J in ['0.5','0','-0.5','-1.0']:
    Ls=[L for L in [17,25,33,49,65] if os.path.exists(f'out2/corner_J{J}_L{L}.bin')]
    if len(Ls)<3: continue
    R=[corner_ratio(load(f'out2/corner_J{J}_L{L}.bin',L),L) for L in Ls]
    m=[i for i,L in enumerate(Ls) if L>=25]
    s,e=jk([R[i] for i in m],np.array([Ls[i] for i in m]))
    out['corner'][J]=dict(L=Ls,ratio=[float(r.mean()) for r in R],err=[float(r.std()/3) for r in R],gamma=-s,gerr=e)
    print('corner J=%s L=%s gamma(L>=25)=%.3f +- %.3f'%(J,Ls,-s,e))
for J in ['-1.0','-0.4']:
    for pos in ['bulk','edge']:
        Ls=[L for L in [16,24,32,48,64] if os.path.exists(f'out2/vac_{pos}_J{J}_L{L}.bin')]
        if not Ls: continue
        P=[]
        for L in Ls:
            site=(L//2-1)+(L//2-1)*L if pos=='bulk' else (L//2 if (L//2)%2==0 else L//2-1)
            P.append(near(load(f'out2/vac_{pos}_J{J}_L{L}.bin',L),L,site))
        Pm=[float(p.mean()) for p in P]; Pe=[float(p.std()/3) for p in P]
        sl=[float(np.log(Pm[i+1]/Pm[i])/np.log(Ls[i+1]/Ls[i])) for i in range(len(Ls)-1)]
        out['vac'][f'{pos}_{J}']=dict(L=Ls,P=Pm,err=Pe,local_slopes=sl)
        print('vac %s J=%s L=%s P=%s +- %s slopes %s'%(pos,J,Ls,np.round(Pm,4),np.round(Pe,4),np.round(sl,3)))
json.dump(out,open('../../data/paper2_mono2.json','w'),indent=1)
