import json, numpy as np
d=json.load(open('../data/paper2_pair_bulk.json'))
print('t     C_inf(r=32)    local xi_eff(r) from D(r)=C(r)-C(r+2)')
res={}
for t in d:
    if t=='0': continue
    a=np.array(d[t]); r=a[:,0]; C=a[:,1]
    D=C[:-2]-C[2:]; rr=r[:-2]
    ok=D>1e-13*C[0]            # stop before round-off floor
    xi=[2/np.log(D[i]/D[i+2]) for i in range(len(D)-2) if ok[i] and ok[i+2] and D[i]/D[i+2]>1]
    xi=np.array(xi)
    res[t]=xi
    print('%-5s %.10f   '%(t,C[-1])+' '.join('%.2f'%x for x in xi[:14:2]))
