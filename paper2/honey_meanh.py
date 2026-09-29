"""Principled corner jump: mean height <h> on faces of balanced armchair flakes.
Crossing bond (a,b) with A on the left: dh = 1 - 3 x_ab  ->  <dh> = 1 - 3 P_ab (boundary bonds: P=0).
<h> integrated on the dual graph (curl-free on average). J_eff(corner) = difference of <h>
averaged over faces at depth d inside the two adjacent edges, away from the corners."""
import numpy as np, sys, collections; sys.path.insert(0,'.'); from honeyflake import flake, D
s3=np.sqrt(3)
def meanh(A,B):
    key=lambda p:(round(p[0],3),round(p[1],3))
    Bi={key(b):j for j,b in enumerate(B)}; Ai={key(a):i for i,a in enumerate(A)}
    K=np.zeros((len(A),len(B))); bonds=[]
    for i,a in enumerate(A):
        for d in D:
            j=Bi.get(key(a+d))
            if j is not None: K[i,j]=1
    G=np.linalg.inv(K)
    adj=collections.defaultdict(list); pos={}
    def dual(a,b,dh):
        m=(a+b)/2; d=b-a; n=np.array([-d[1],d[0]])
        P,Q=m+n*s3/2,m-n*s3/2
        if (Q-P)[0]*(a-P)[1]-(Q-P)[1]*(a-P)[0]<0: P,Q=Q,P      # A on the left going P->Q
        kp,kq=key(P),key(Q); pos[kp]=P; pos[kq]=Q
        adj[kp].append((kq,dh)); adj[kq].append((kp,-dh))
    for i,a in enumerate(A):
        for d in D:
            b=a+d; j=Bi.get(key(b))
            dual(a,b, 1-3*K[i,j]*G[j,i] if j is not None else 1.0)
    for j,b in enumerate(B):                     # boundary bonds seen from B (A missing)
        for d in D:
            a=b-d
            if key(a) not in Ai: dual(a,b,1.0)
    root=min(pos,key=lambda k:np.linalg.norm(pos[k]-np.array([0,-1.0]))); h={root:0.0}; q=[root]
    while q:
        u=q.pop()
        for v,dh in adj[u]:
            if v not in h: h[v]=h[u]+dh; q.append(v)
    P=np.array([pos[k] for k in h]); H=np.array([h[k] for k in h])
    # consistency: max curl error
    err=max(abs(h[v]-h[u]-dh) for u in adj for v,dh in adj[u] if u in h and v in h)
    return P,H,err
def edge_levels(P,H,R,depths,c0=np.array([0,-1.0]),frac=0.5):
    out=[]
    for k in range(6):
        n=np.array([np.cos(np.radians(60*k)),np.sin(np.radians(60*k))]); t=np.array([-n[1],n[0]])
        X=P-c0; dep=R-X@n; along=X@t; half=R*np.tan(np.pi/6)
        row=[]
        for d0,d1 in depths:
            m=(dep>=d0)&(dep<d1)&(np.abs(along)<frac*half); row.append(H[m].mean() if m.any() else np.nan)
        out.append(row)
    return np.array(out)      # [edge k][depth bin]
if __name__=='__main__':
    depths=[(-1.5,0.2),(0.2,2),(2,4),(4,6),(6,8)]
    for R in [int(x) for x in sys.argv[1:]]:
        A,B=flake(R,'armchair')
        if len(A)!=len(B): print(R,'unbalanced'); continue
        P,H,err=meanh(A,B); L=edge_levels(P,H,R,depths)
        print('R=%d sites=%d curl-err=%.1e'%(R,2*len(A),err))
        for di,(d0,d1) in enumerate(depths):
            lev=L[:,di]; J=[lev[k]-lev[k-1] for k in range(6)]
            print('   depth [%4.1f,%4.1f): edge <h> %s | J at corners 30..330: %s'%(d0,d1,' '.join('%+.3f'%x for x in lev),' '.join('%+.3f'%x for x in J)))
