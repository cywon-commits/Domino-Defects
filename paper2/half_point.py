import sys, numpy as np; sys.path.insert(0,'.'); from pftri import logZ, logZ_multi
L=48; t=0.5; lz=logZ(L,L,t); c0=L//2; q=np.sqrt(3)-2
print('horizontal formula (1-(sqrt3-2)^r)/6')
for name,step in [('vertical (dr,dc)=(1,0)',(1,0)),('diagonal NE (1,1)',(1,1)),('anti-diagonal (1,-1)',(1,-1)),('knight (1,2)',(1,2))]:
    s=[]
    for k in range(1,9):
        dr,dc=step[0]*k,step[1]*k
        u=(c0-dr//2,c0-dc//2); v=(u[0]+dr,u[1]+dc)
        C=np.exp(logZ_multi(L,L,t,[u,v])-lz); s.append(C)
    s=np.array(s); d=s-1/6
    print('%-24s C-1/6: %s   ratios: %s'%(name,' '.join('%+.3e'%x for x in d[:5]),' '.join('%.5f'%(d[i+1]/d[i]) for i in range(5))))
