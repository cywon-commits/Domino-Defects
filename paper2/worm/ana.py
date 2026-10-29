import numpy as np, sys, os, json
W,H=64,128
def row_corr(fn,row):
    h=np.fromfile(fn).reshape(H,W)[row]; x=np.arange(1,W,2); c=h[x]   # opposite colour to tail at x=0
    l=(W/np.pi)*np.sin(np.pi*x/W); return x,l,c
def slope(l,c,lo,hi):
    m=(l>=lo)&(l<=hi)&(c>0); p,cov=np.polyfit(np.log(l[m]),np.log(c[m]),1,cov=True); return p[0],np.sqrt(cov[0,0])
res={}
for J in sys.argv[1:]:
    fb='out/J%s_bulk.bin'%J; fe='out/J%s_edge.bin'%J
    if not (os.path.exists(fb) and os.path.exists(fe)): continue
    xb,lb,cb=row_corr(fb,64); xe,le,ce=row_corr(fe,0)
    out=[]
    for lo,hi in [(4,20),(6,20),(8,20.4)]:
        sb,eb=slope(lb,cb,lo,hi); se,ee=slope(le,ce,lo,hi); out.append((sb,se,se/sb))
    res[J]=out
    print('beta*v=%5s  bulk slope (=-g)   edge slope (=-2g pred)   ratio   [chord windows 4-20, 6-20, 8-20]'%J)
    for (sb,se,r),w in zip(out,['4-20','6-20','8-20']): print('        %s: %+.3f   %+.3f   %.3f'%(w,sb,se,r))
json.dump(res,open('ana.json','w'))
