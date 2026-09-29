"""delta from a long balanced armchair strip between edges with normals 0 (class +) and 180 (class -)."""
import numpy as np, sys, json; sys.path.insert(0,'.')
import honey_defect as hd; from honey_corner import boundary; from honey_meanh import meanh
s3=np.sqrt(3); D=hd.D; n0=np.array([1.0,0]); n60=np.array([0.5,s3/2])
def strip(W,L,c0):
    a1=np.array([s3,0]); a2=np.array([s3/2,1.5]); M=int(1.4*(W+L))+6
    A=np.array([m*a1+n*a2 for m in range(-M,M+1) for n in range(-M,M+1)]); B=A+D[0]
    ins=lambda P:(np.abs((P-c0)@n0)<=W)&(np.abs((P-c0)@n60)<=L); A=A[ins(A)]; B=B[ins(B)]
    key=lambda p:(round(p[0],4),round(p[1],4))
    while True:
        Bs={key(b) for b in B}; As={key(a) for a in A}
        dA=np.array([sum(key(a+d) in Bs for d in D) for a in A]); dB=np.array([sum(key(b-d) in As for d in D) for b in B])
        if (dA>=2).all() and (dB>=2).all(): break
        A=A[dA>=2]; B=B[dB>=2]
    return A,B
out=[]
for W in [6,8,10,12]:
    L=3.5*W
    for cname,c0 in [('A',np.zeros(2)),('face',np.array([0,-1.0])),('B',np.array([0,1.0]))]:
        A,B=strip(W,L,c0)
        if len(A)!=len(B): continue
        X,Hb,tot=boundary(A,B); Pb=X-c0; xb=Pb@n0; yb=Pb@n60
        mid=np.abs(yb)<0.3*L
        e0=mid&(xb>0.8*xb.max()); e180=mid&(xb<0.8*xb.min())
        m0,m180=np.median(Hb[e0]),np.median(Hb[e180]); x0,x180=xb[e0].mean(),xb[e180].mean()
        P,H,err=meanh(A,B); Pc=P-c0; x=Pc@n0; y=Pc@n60
        sel=(np.abs(y)<0.3*L)&(x<x0-3)&(x>x180+3)
        beta,alpha=np.polyfit(x[sel],H[sel],1)
        b0,b180=alpha+beta*x0,alpha+beta*x180
        rec=dict(W=W,centre=cname,nA=len(A),m0=m0,m180=m180,delta0=b0-m0,delta180=b180-m180,slope=beta)
        out.append(rec)
        print('W=%2d centre=%-4s nA=%4d  midpoint levels m0=%+.2f m180=%+.2f  -> delta(class+,n=0)=%+.4f  delta(class-,n=180)=%+.4f  (slope %.4f, curl %.0e)'%(W,cname,len(A),m0,m180,b0-m0,b180-m180,beta,err),flush=True)
json.dump(out,open('../data/paper2_honey_strip.json','w'),indent=1)
