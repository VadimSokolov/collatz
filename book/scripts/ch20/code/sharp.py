"""Sharpness of the cycle inequality: rational cycles of Christoffel words.
For s/r on the ladder, the lower Christoffel word K_t = floor(t s / r) defines a
rational cycle. We compute its minimum M, the ratio M/F(r) with F = r/(3 Lambda),
and the ratio 3 M Lambda / B_r, where B_r = sum_{t<r} 2^{-{t theta}}."""
import numpy as np, mpmath as mp, math
from ladder import rungs, TH
mp.mp.dps = 60
L2 = math.log(2)
th = float(TH)
def christoffel_cycle(r, s):
    t = np.arange(r, dtype=np.int64)
    K = (t * s) // r                       # exact integer floors
    # W_t = K_t - t theta, computed stably: t*theta split into int and frac
    ft = np.array([float(mp.frac(int(x)*TH)) for x in t]) if r < 3000 else None
    if ft is None:
        # t*theta = t + t*(theta-1); use high/low split of theta-1
        hi = float(TH - 1); lo = float((TH - 1) - mp.mpf(hi))
        prod = t.astype(np.float64)*hi
        ip = np.floor(prod); fr = prod - ip + t.astype(np.float64)*lo
        ft = fr - np.floor(fr)
        whole = t + ip.astype(np.int64) + np.floor(prod - ip + t*lo).astype(np.int64)
    else:
        whole = np.array([int(mp.floor(int(x)*TH)) for x in t])
    W = (K - whole).astype(np.float64) - ft   # K_t - t theta
    a = np.exp2(W)                          # a_t = 2^{K_t}/3^t
    d = s - r*TH                            # delta
    Lam = float(d*mp.log(2))
    x0 = a.sum()/(3*float(mp.expm1(d*mp.log(2))))   # x_0 = sum a_j / (3(2^delta - 1))
    cum = np.concatenate(([0.0], np.cumsum(a)[:-1]))
    x = (x0 + cum/3)/a                      # x_t
    M = x.min()
    B = np.exp2(-ft).sum()
    return M, Lam, B, x.max()
if __name__ == "__main__":
    print(f"{'r':>10} {'s':>10} {'M/F':>10} {'3M Lam/B_r':>12} {'xmax/F':>9} {'B_r/r':>10}")
    for r, s, d in rungs:
        if r > 2e7: break
        M, Lam, B, X = christoffel_cycle(r, s)
        F = r/(3*Lam)
        print(f"{r:>10d} {s:>10d} {M/F:10.6f} {3*M*Lam/B:12.6f} {X/F:9.5f} {B/r:10.6f}")
    print("1/(2 ln 2) =", 1/(2*L2))
