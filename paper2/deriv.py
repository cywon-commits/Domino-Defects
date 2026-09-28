"""Small-t expansion of ln R(N,t) - ln R(N,0) = a1(N) t + a2(N) t^2 : finite differences."""
import sys, json, time, numpy as np; sys.path.insert(0,'.'); from pftri import logZ
h=2e-3; out=[]; T0=time.time()
for N in [13,21,29,37,45,53,61,69,81]:
    m=(N-1)//2
    def lr(t,c): return logZ(N,N,t,c)-logZ(N,N,t,(m,m))
    row=[N]
    for c in [(0,0),(0,N-1)]:
        f0,fp,fp2=lr(0,c),lr(h,c),lr(2*h,c)
        a1=(-3*f0+4*fp-fp2)/(2*h); a2=(f0-2*fp+fp2)/(2*h*h)
        row+= [a1,a2]
    out.append(row); print('N=%d  A: a1=%+.5f a2=%+.3f | B: a1=%+.5f a2=%+.3f  (%.0fs)'%(N,*row[1:],time.time()-T0),flush=True)
json.dump(out,open('../data/paper2_deriv.json','w'))
a=np.array(out); N=a[:,0]
for k,name in [(2,'A a2'),(4,'B a2')]:
    p=np.polyfit(np.log(N[3:]),np.log(np.abs(a[3:,k])),1); print(name,'~ N^%.3f'%p[0])
