"""Boundary height of honeycomb flakes (lozenge Thurston units: step +1 if the interior vertex
next to the boundary edge is A, -1 if B), average height per straight edge, corner jumps J_c."""
import numpy as np, sys; sys.path.insert(0,'.'); from honeyflake import flake, D
def boundary(A,B):
    key=lambda p:(round(p[0],3),round(p[1],3))
    inA={key(a) for a in A}; inB={key(b) for b in B}; edges={}
    def add(v,u,sgn):
        m=(v+u)/2; d=u-v; n=np.array([-d[1],d[0]])/np.linalg.norm(d)
        P,Q=m+n*np.sqrt(3)/2, m-n*np.sqrt(3)/2
        if (Q-P)[0]*(v-P)[1]-(Q-P)[1]*(v-P)[0] < 0: P,Q=Q,P      # interior v on the left
        edges[key(P)]=(key(Q),sgn,(P+Q)/2)
    for a in A:
        for d in D:
            if key(a+d) not in inB: add(a,a+d,+1)
    for b in B:
        for d in D:
            if key(b-d) not in inA: add(b,b-d,-1)
    start=next(iter(edges)); cur=start; h=0; H=[]; X=[]
    for _ in range(len(edges)):
        nxt,s,mid=edges[cur]; H.append(h+s/2); X.append(mid); h+=s; cur=nxt
        if cur==start: break
    assert cur==start and len(H)==len(edges), 'boundary not a single cycle'
    return np.array(X),np.array(H),h
if __name__=="__main__":
  for kind,R in [('armchair',12),('armchair',20),('armchair',28),('zigzag',12)]:
    A,B=flake(R,kind); X,H,tot=boundary(A,B); c=np.array([0,-1.0])
    ang=np.degrees(np.arctan2(X[:,1]-c[1],X[:,0]-c[0]))%360
    corners=np.arange(30,360,60) if kind=='armchair' else np.arange(0,360,60)
    # assign each boundary step to the edge (sector between consecutive corners), skip 12 deg near corners
    avg=[]
    for k in range(6):
        lo,hi=corners[k],corners[(k+1)%6]+(360 if k==5 else 0)
        a2=np.where(ang<corners[0],ang+360,ang)
        m=(a2>lo+12)&(a2<hi-12); avg.append(H[m].mean())
    J=[avg[k]-avg[k-1] for k in range(6)]
    print('%-8s R=%d  boundary steps=%d  net winding=%d  edge averages: %s   corner jumps J (at %s deg): %s'%(kind,R,len(H),tot,
          ' '.join('%+.3f'%v for v in avg), ','.join(str(x) for x in corners), ' '.join('%+.3f'%v for v in J)))
