import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
# Where is the variance deficit of the tail parity sum?  Split each orbit of n in a dyadic
# block at the first time it falls below Y = 2^y: middle part (after the Terras horizon k,
# before falling below Y) and end part (from below Y to 1, read from the census table).
import numpy as np, sys
tau = np.load(BOOK + '/scripts/data/ch3_tau_small.npy').astype(np.int64)  # n < 2^26
L3 = np.log(3.0); S2 = (L3 / 2) ** 2
# odd counts for n < 2^20 recomputed directly
def odd_count(x):
    a = 0
    while x > 1:
        if x & 1: x = (3 * x + 1) >> 1; a += 1
        else: x >>= 1
    return a
rng = np.random.default_rng(5)
for k in (20, 24, 28):
    for y in (6, 10, 14):
        Y = 1 << y
        ns = rng.integers(1 << k, 1 << (k + 1), 200000)
        Mh = []; Mm = []; Me = []; Tm = []; Te = []
        for n in ns:
            x = int(n); s = 0; ah = 0; am = 0
            while s < k:
                if x & 1: x = (3 * x + 1) >> 1; ah += 1
                else: x >>= 1
                s += 1
            sm = 0
            while x >= Y:
                if x & 1: x = (3 * x + 1) >> 1; am += 1
                else: x >>= 1
                sm += 1
            te = int(tau[x]); ae = odd_count(x)
            Mh.append(L3 / 2 * (2 * ah - k)); Mm.append(L3 / 2 * (2 * am - sm)); Me.append(L3 / 2 * (2 * ae - te))
            Tm.append(sm); Te.append(te)
        Mh, Mm, Me, Tm, Te = map(np.array, (Mh, Mm, Me, Tm, Te))
        print(f"k={k} Y=2^{y}: Var Mh/(s2 k)={Mh.var()/(S2*k):.4f}  Var Mmid/(s2 E steps)={Mm.var()/(S2*Tm.mean()):.4f} "
              f"(E steps {Tm.mean():.1f})  Var Mend={Me.var():.3f} vs s2 E steps={S2*Te.mean():.3f} (E {Te.mean():.1f}) "
              f"deficit end={Me.var()-S2*Te.mean():.3f}  cov(mid,end)={np.cov(Mm,Me)[0,1]:.3f} cov(h,mid)={np.cov(Mh,Mm)[0,1]:.3f}", flush=True)
