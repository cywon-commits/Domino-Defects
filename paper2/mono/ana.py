import numpy as np,glob,json,sys
def ratio(fn,L,r=2):
    h=np.fromfile(fn).reshape(L,L); h/=h.sum()
    x,y=np.meshgrid(range(L),range(L),indexing='ij')
    black=((x+y)%2==0)
    c=0.;n=0
    for cx,cy in [(0,0),(0,L-1),(L-1,0),(L-1,L-1)]:
        m=black&(abs(x-cx)+abs(y-cy)<=r); c+=h[m].sum(); n+=m.sum()
    mc=black&(abs(x-L//2)+abs(y-L//2)<=r)
    return (c/n)/(h[mc].sum()/mc.sum())
out={}
for J in ['0.5','0','-0.5','-1.0']:
    Ls=[];R=[]
    for L in [17,25,33,49,65]:
        fn=f'out/J{J}_L{L}.bin'
        try: R.append(ratio(fn,L)); Ls.append(L)
        except Exception: pass
    if len(Ls)>=2:
        sl=np.polyfit(np.log(Ls),np.log(R),1)[0]
        out[J]=dict(L=Ls,ratio=R,gamma=sl)
        print(J,Ls,np.round(R,4),'gamma',round(sl,3))
json.dump(out,open('corner.json','w'),indent=1)
