import mpmath as mp, sys; sys.path.insert(0,'.'); from bb import coeffs, K2
co=coeffs('0.5',31,dps=50); mp.mp.dps=40; q=mp.sqrt(3)-2
worst=max(abs(K2(co,n)-(1-q**(n+1))/6) for n in range(1,31)); print('t=1/2: max |K2(n)-(1-q^(n+1))/6| over n=1..30 =',mp.nstr(worst,3))
# t=0.3 control: does D_n=4K2^2 satisfy any linear recurrence of order <=4? test via Hankel determinants of D_n
co=coeffs('0.3',26,dps=50); mp.mp.dps=40
D=[4*K2(co,n)**2 for n in range(1,26)]
for r in [2,3,4,5]:
    H=mp.matrix(r+1,r+1)
    for i in range(r+1):
        for j in range(r+1): H[i,j]=D[10+i+j]
    print('t=0.3: Hankel det of order %d (zero iff a recurrence of order %d holds):'%(r+1,r), mp.nstr(mp.det(H),3))
co=coeffs('0.5',26,dps=50); mp.mp.dps=40
D=[4*K2(co,n)**2 for n in range(1,26)]
for r in [2,3,4]:
    H=mp.matrix(r+1,r+1)
    for i in range(r+1):
        for j in range(r+1): H[i,j]=D[10+i+j]
    print('t=1/2: Hankel det of order %d:'%(r+1), mp.nstr(mp.det(H),3))
