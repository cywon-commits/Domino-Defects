import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
import sys, math, time
sys.argv=['x']; exec(open('lgv2.py').read().split("if __name__")[0])
import mpmath as mp
H=lambda x: -(x*mp.log(x)+(1-x)*mp.log(1-x)) if 0<x<1 else mp.mpf(0)
def predict(al,be,ga,kap,lam,cc):
    fr=lambda l: (al+l)*H(l/(al+l)) + (ga+be-l)*H((be-l)/(ga+be-l))
    fb=lambda k: (be+k)*H(k/(be+k)) + (ga+al-k)*H((al-k)/(ga+al-k))
    Lr=lambda l: fr(l)-fr(0); Lb=lambda k: fb(k)-fb(0)
    g=lambda t: t*(Lb(kap/t)+Lr(lam/t)-cc)
    tmin=max(kap/al,lam/be); ts=[tmin+(1-tmin)*i/3000 for i in range(1,3001)]
    b=max(ts,key=g); return float(g(1.0)), float(g(b)), b
def dlogZc(al,be,ga,N=600):
    a,b,c=int(al*N),int(be*N),int(ga*N)
    from functools import lru_cache
    def lG(n):
        s=mp.mpf(0); f=mp.mpf(0)
        for k in range(1,n): f+=mp.log(k); s+=f
        return s
    lZ=lambda a,b,c: lG(a)+lG(b)+lG(c)+lG(a+b+c)-lG(a+b)-lG(b+c)-lG(c+a)
    return float((lZ(a,b,c+1)-lZ(a,b,c))/N)
cases=[(1,1,1),(1,1.5,0.75)]
for al,be,ga in cases:
    cc=dlogZc(al,be,ga)
    print('hexagon',(al,be,ga),'layer cost d_c logZ/n = %.5f'%cc)
    for kap,lam in ((0.1,0.1),(0.1,0.3),(0.2,0.2),(0.3,0.3),(0.45,0.5)):
        p1,pm,tb=predict(al,be,ga,kap,lam,cc)
        row=[]
        for n in (60,120,240):
            a,b,c=int(al*n),int(be*n),int(ga*n)
            k=int(kap*n)+1; l=int(lam*n)+1
            t=time.time(); w=lgvW(a,b,c,k,l)/lgvZ(a,b,c)
            row.append(float(mp.log(w)/n))
        print('  (%.2f,%.2f) exact n=60,120,240: %s | rule tau=1 %.4f, rule max %.4f (tau*=%.2f)'%(kap,lam,[round(r,4) for r in row],p1,pm,tb),flush=True)
