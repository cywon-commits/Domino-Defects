"""Model-free test: Y_p(r) = (C(r)-Cinf) e^{2 kappa r} for each parity p; local exponent
s_p(r) = dlnY/dlnr should tend to a constant (power-law prefactor) iff the rate is exactly 2kappa.
Also report the residual rate delta_p(r) = -d ln Y/dr."""
import json,sys,numpy as np; sys.path.insert(0,'.'); from symbol_xi import kappa
d=json.load(open('../data/paper2_pair_big.json'))
for t in sys.argv[1:]:
    C=np.array(d[t]['C']); r=C[:,0].astype(int); C=C[:,1]; k=kappa(float(t),2001)
    # Cinf: extrapolate each parity's last points geometrically (rate 4k per step of 2)
    q=np.exp(-4*k); Cinf=np.mean([C[-1]-(C[-3]-C[-1])*q/(1-q), C[-2]-(C[-4]-C[-2])*q/(1-q)])
    P=C-Cinf; line='t=%s 2k=%.4f Cinf=%.13f\n'%(t,2*k,Cinf)
    for par in [1,0]:
        m=(r%2==par)&(np.abs(P)>2e-12)
        rr=r[m]; Y=P[m]*np.exp(2*k*rr)
        s=np.diff(np.log(np.abs(Y)))/np.diff(np.log(rr)); dl=-np.diff(np.log(np.abs(Y)))/np.diff(rr)
        sign=''.join('+' if y>0 else '-' for y in Y[:6])
        pick=[i for i in range(0,len(s),max(1,len(s)//8))]
        line+='   %s r: sign %s | local exponent s(r): %s\n'%('odd ' if par else 'even',sign,' '.join('%d:%.2f'%((rr[i]+rr[i+1])/2,s[i]) for i in pick))
        line+='        residual rate delta(r)=%s\n'%(' '.join('%.4f'%dl[i] for i in pick))
    print(line)
