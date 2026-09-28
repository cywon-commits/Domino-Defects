import sys, json, time, numpy as np; sys.path.insert(0,'.'); from pftri import logZ
a=np.array(json.load(open('../data/paper2_deriv.json'))); N=a[:,0]; m=N>=29
preds={}
for k,n in [(2,'A'),(4,'B')]:
    for lab,f in [('N2+N',lambda x:np.c_[x**2,x,np.ones_like(x)]),('N2+NlnN',lambda x:np.c_[x**2,x*np.log(x),np.ones_like(x)]),('pure N^p',None)]:
        if f is None:
            p=np.polyfit(np.log(N[m]),np.log(a[m,k]),1); preds[(n,lab)]=np.exp(np.polyval(p,np.log(101)))
        else:
            c,*_=np.linalg.lstsq(f(N[m]),a[m,k],rcond=None); preds[(n,lab)]=(f(np.array([101.]))@c)[0]
T0=time.time(); Nn=101; mm=50; h=2e-3
def lr(t,c): return logZ(Nn,Nn,t,c)-logZ(Nn,Nn,t,(mm,mm))
for n,c in [('A',(0,0)),('B',(0,Nn-1))]:
    f0,fp,fp2=lr(0,c),lr(h,c),lr(2*h,c); a2=(f0-2*fp+fp2)/(2*h*h)
    print(n,'measured a2(101)=%.3f'%a2,' predictions:',', '.join('%s %.3f'%(lab,preds[(n,lab)]) for lab in ['N2+N','N2+NlnN','pure N^p']),' (%.0fs)'%(time.time()-T0),flush=True)
