"""Exact count of admissible words of length r (K_r = s) by the number n of indices
t < r with W_t = K_t - t theta > -b.  Forward DP over (K_t, n) in two layers, counts
kept as floating point with a running log scale.  Output: log10 of the number of
admissible words with n >= f r, for a grid of f, and the total."""
import numpy as np, mpmath as mp, math, sys
mp.mp.dps = 40
TH = mp.log(3, 2)
def count(r, s, b):
    d = s - r*TH
    cap = [int(mp.floor(t*TH + d)) for t in range(r)]
    low = [int(mp.floor(t*TH - b)) + 1 for t in range(r)]    # W_t > -b  <=>  K_t >= low[t]
    old = np.zeros((s+1, r+1)); old[0, 1 if 0 >= low[0] else 0] = 1.0
    L = np.zeros(r+1)                                          # per column log10 scale
    Ks = np.arange(s+1)
    for t in range(1, r):
        C = np.cumsum(old, axis=0)
        S = np.zeros_like(old); S[1:] = C[:-1]                 # sum over K < K'
        inb = Ks >= low[t]
        A = S.copy(); A[inb] = 0.0                             # stay in column n
        B = np.zeros_like(S); B[inb, 1:] = S[inb, :-1]         # move to column n+1
        LB = np.concatenate(([-np.inf], L[:-1]))
        M = np.maximum(L, LB)
        with np.errstate(invalid="ignore"):
            new = A*10.0**(np.where(np.isfinite(M), L - M, 0.0)) + B*10.0**(np.where(np.isfinite(LB), LB - M, -np.inf))
        new[cap[t]+1:] = 0.0; new[:t] = 0.0
        m = new.max(axis=0)
        with np.errstate(divide="ignore"):
            L = np.where(m > 0, M + np.log10(np.maximum(m, 1e-300)), -np.inf)
        new = np.where(m > 0, new/np.maximum(m, 1e-300), 0.0)
        old = new
    N = old.sum(axis=0)
    with np.errstate(divide="ignore"):
        logN = np.where(N > 0, np.log10(np.maximum(N, 1e-300)) + L, -np.inf)
    m = logN.max()
    return logN, m + math.log10(np.sum(10.0**(logN - m)))
if __name__ == "__main__":
    rungs = {41: 65, 306: 485, 971: 1539, 2301: 3647}
    if len(sys.argv) > 1: rungs = {int(sys.argv[1]): {41: 65, 306: 485, 971: 1539, 2301: 3647}[int(sys.argv[1])]}
    fs = [0.084, 0.2, 0.3, 0.39, 0.476, 0.55, 0.65, 0.75, 0.9]
    for b in (1, 3):
        print(f"\nband b = {b}: log10 of the number of admissible words with at least f r indices in band, divided by r")
        print("      r  total/r " + " ".join(f"f={f:<5}" for f in fs))
        for r, s in rungs.items():
            logN, logtot = count(r, s, b)
            tail = []
            for f in fs:
                n0 = math.ceil(f*r)
                v = logN[n0:]; m = v.max()
                tail.append((m + math.log10(np.sum(10.0**(v - m))))/r if np.isfinite(m) else -np.inf)
            np.save(f"../data/bandcount_b{b}_r{r}.npy", logN)
            print(f"{r:7d}  {logtot/r:7.4f} " + " ".join(f"{x:7.4f}" for x in tail), flush=True)
