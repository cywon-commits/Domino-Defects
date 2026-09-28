"""Task 1 step 3: corner/center ratio R = Z(board - corner)/Z(board - center), odd N x N."""
import sys, json, time, numpy as np
sys.path.insert(0, '.'); from pftri import logZ
ts = [0, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]
Ns = [11, 15, 21, 31, 41, 51, 61]
res = {}
t0 = time.time()
for t in ts:
    row = []
    for N in Ns:
        m = (N-1)//2
        lc = logZ(N, N, t, (0, 0)); lm = logZ(N, N, t, (m, m))
        # other corner: (0,N-1) is inequivalent to (0,0) because of the diagonal direction
        lo = logZ(N, N, t, (0, N-1))
        row.append((N, float(np.exp(lc-lm)), float(np.exp(lo-lm))))
    res[str(t)] = row
    print(t, ['%d:%.4g/%.4g' % r for r in row], '%.0fs' % (time.time()-t0), flush=True)
json.dump(res, open('../data/paper2_ratio.json', 'w'), indent=1)
