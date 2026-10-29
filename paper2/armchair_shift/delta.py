import numpy as np, sys, collections, json; sys.path.insert(0,'.'); from armstrip import build, probs, D, s3
def heights(A,B,bonds,P):
    key=lambda p:(round(p[0],3),round(p[1]%3,3))
    Ak={key(a) for a in A}; Bk={key(b) for b in B}
    Pb={}
    for (i,j,s,t),p in zip(bonds,P): Pb[(key(A[i]),t)]=p
    adj=collections.defaultdict(list); pos={}; bnd=[]
    def dual(a,d,p,present_both):
        b=a+d; m=(a+b)/2; n=np.array([-d[1],d[0]])
        Pp,Q=m+n*s3/2,m-n*s3/2
        if (Q-Pp)[0]*(a-Pp)[1]-(Q-Pp)[1]*(a-Pp)[0]<0: Pp,Q=Q,Pp      # A on the left going P->Q
        kp,kq=key(Pp),key(Q); pos[kp]=Pp; pos[kq]=Q
        dh=1-3*p; adj[kp].append((kq,dh,Q[1]-Pp[1])); adj[kq].append((kp,-dh,Pp[1]-Q[1]))
        if not present_both: bnd.append((kp,kq,m))
    for a in A:
        for t,d in enumerate(D):
            if key(a+d) in Bk: dual(a,d,Pb[(key(a),t)],True)
            else: dual(a,d,0.0,False)
    for b in B:
        for t,d in enumerate(D):
            if key(b-d) not in Ak: dual(b-d,d,0.0,False)
    root=next(iter(pos)); h={root:0.0}; yy={root:pos[root][1]}; q=[root]; err=0
    while q:
        u=q.pop()
        for v,dh,dy in adj[u]:
            if v not in h: h[v]=h[u]+dh; q.append(v)
            else: err=max(err,abs(h[v]-h[u]-dh))
    return h,pos,bnd,err
def delta(R,nk=256,dcut=6,pre=None):
    A,B,bonds=build(R) if pre is None else pre[:3]
    P=probs(A,B,bonds,nk)[0] if pre is None else pre[3]
    h,pos,bnd,err=heights(A,B,bonds,P)
    Ak={(round(a[0],3),round(a[1]%3,3)) for a in A}; Bk={(round(b[0],3),round(b[1]%3,3)) for b in B}
    # interior faces: all six vertices present
    X=[];Hh=[]
    for k,p in pos.items():
        verts=[p+np.array([0,1.0]),p+np.array([s3/2,0.5]),p+np.array([s3/2,-0.5]),p+np.array([0,-1.0]),p+np.array([-s3/2,-0.5]),p+np.array([-s3/2,0.5])]
        kk=[(round(v[0],3),round(v[1]%3,3)) for v in verts]
        if all(x in Ak or x in Bk for x in kk): X.append(p[0]); Hh.append(h[k])
    X=np.array(X); Hh=np.array(Hh); xmax=X.max(); xmin=X.min()
    sel=(X<xmax-dcut)&(X>xmin+dcut); sl,ic=np.polyfit(X[sel],Hh[sel],1)
    mR=[(h[a]+h[b])/2 for a,b,m in bnd if m[0]>0]; xR=np.mean([m[0] for a,b,m in bnd if m[0]>0])
    mL=[(h[a]+h[b])/2 for a,b,m in bnd if m[0]<0]; xL=np.mean([m[0] for a,b,m in bnd if m[0]<0])
    dR=ic+sl*xR-np.mean(mR); dL=ic+sl*xL-np.mean(mL)
    return dict(R=R,nA=len(A),curl=err,spreadR=float(np.ptp(mR)),spreadL=float(np.ptp(mL)),slope=float(sl),dplus=float(dR),dminus=float(dL),diff=float(dR-dL))
if __name__=='__main__':
    out=[]
    for R in [float(x) for x in sys.argv[1:]]:
        r=delta(R); out.append(r)
        print('R=%5g nA=%4d curl=%.0e  m-spread R/L=%.2g/%.2g  slope=%.2e  delta(+)=%+.6f delta(-)=%+.6f  diff=%+.6f'%(r['R'],r['nA'],r['curl'],r['spreadR'],r['spreadL'],r['slope'],r['dplus'],r['dminus'],r['diff']),flush=True)
    json.dump(out,open('delta_%s.json'%sys.argv[1],'w'))
