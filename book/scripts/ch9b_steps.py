"""Step counts for the chapter on positive density from the backward tree (adapted from N. Polson's collatz_steps.py).
For 2 <= n <= 3*10^6: sigma0(n), the number of steps of the standard map T_0 to reach 1, and A(n), the number of odd steps
(the number of steps of the Syracuse map). Per range of n: mean and standard deviation of sigma0/log n, mean of A/log n, the
share meeting Lau's bounds sigma0 <= floor(10.46 log n) and A <= floor(3.4881 log n), and how many n meet exactly one of them.
Writes data/ch9b_steps.json, with the histograms of sigma0/log n and A/log n on [10^6, 3*10^6) for ch9b_figures.py."""
import os, json, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, 'data'); os.makedirs(DATA, exist_ok=True)
N = 3_000_000
tot = np.zeros(N + 1, dtype=np.int32); acc = np.zeros(N + 1, dtype=np.int32)
for n in range(2, N + 1):
    m, t, a = n, 0, 0
    while m >= n:                      # run until the orbit drops below n, then reuse the stored counts
        if m % 2 == 0: m //= 2
        else: m = 3 * m + 1; a += 1
        t += 1
    tot[n] = t + tot[m]; acc[n] = a + acc[m]
ns = np.arange(2, N + 1); L = np.log(ns); T = tot[2:]; A = acc[2:]
res = dict(N=N, ranges={})
for lo, hi in [(10**4, 10**5), (10**5, 10**6), (10**6, N)]:
    s = (ns >= lo) & (ns < hi)
    bt = T[s] <= np.floor(10.46 * L[s]); ba = A[s] <= np.floor(3.4881 * L[s])
    res['ranges'][f'{lo}-{hi}'] = dict(count=int(s.sum()), meanT=float((T[s] / L[s]).mean()), sdT=float((T[s] / L[s]).std()),
        meanA=float((A[s] / L[s]).mean()), sdA=float((A[s] / L[s]).std()), fracT=float(bt.mean()), fracA=float(ba.mean()),
        fracBoth=float((bt & ba).mean()), onlyT=int((bt & ~ba).sum()), onlyA=int((ba & ~bt).sum()),
        logmid=math.log((lo + hi) / 2))
s = ns >= 10**6
hT, eT = np.histogram(T[s] / L[s], bins=120, density=True); hA, eA = np.histogram(A[s] / L[s], bins=120, density=True)
res['hist'] = dict(T=dict(dens=hT.tolist(), edges=eT.tolist()), A=dict(dens=hA.tolist(), edges=eA.tolist()))
json.dump(res, open(os.path.join(DATA, 'ch9b_steps.json'), 'w'))
print(json.dumps(res['ranges'], indent=1))
