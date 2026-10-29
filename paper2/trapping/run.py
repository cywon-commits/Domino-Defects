import sys, json, time, numpy as np; sys.path.insert(0,'.'); from trap import Board
def cross(r,c): return [(r,c),(r-1,c),(r+1,c),(r,c-1),(r,c+1)]
def config(case,L):
    m=L//2
    if case=='bulk_q1':  ctr=(m,m);   fixed=[(m,m)]                       # black (m+m even)
    if case=='edge_q1':  ctr=(0,m);   fixed=[(0,m)]                       # black on row 0 (m even)
    if case=='bulk_q3':  ctr=(m,m+1); fixed=cross(m,m+1)+[(L//4,L//4+1),(3*L//4,3*L//4+1)]
    if case=='edge_q3':  ctr=(1,m); fixed=cross(1,m)+[(m,L//4+1),(m,3*L//4+1)]
    return ctr,fixed
out=json.load(open('res.json')) if len(sys.argv)>2 else {}
for L in [int(x) for x in sys.argv[1].split(',')]:
    B=Board(L)
    for case in ['bulk_q1','edge_q1','bulk_q3','edge_q3']:
        ctr,fixed=config(case,L); T0=time.time()
        assert all(0<=r<L and 0<=c<L for r,c in fixed)
        q=sum(1 if (r+c)%2==0 else -1 for r,c in fixed)
        W=[(r,c) for r in range(L) for c in range(L) if (r+c)%2 and (r,c) not in fixed]
        lz=np.array([B.logZ(fixed+[v]) for v in W]); p=np.exp(lz-lz.max()); p/=p.sum()
        d=np.array([np.hypot(r-ctr[0],c-ctr[1]) for r,c in W])
        rec=dict(q_total_removed=q,P_near3=float(p[d<=3.01].sum()),P_near5=float(p[d<=5.01].sum()),
                 prof=[[float(x) for x in d[d<=12]],[float(x) for x in p[d<=12]]])
        out.setdefault(case,{})[str(L)]=rec; json.dump(out,open('res.json','w'))
        print('L=%d %-8s Q_removed=%+d  P(d<=3)=%.4g  P(d<=5)=%.4g  (%.0fs)'%(L,case,q,rec['P_near3'],rec['P_near5'],time.time()-T0),flush=True)
