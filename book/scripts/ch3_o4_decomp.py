"""Chapter 3, problem O4: exact census decomposition of total stopping times.

For every n < 2^LOGN we memoise the total stopping time tau(n) of the
accelerated map T and the number a(n) of odd steps on the way to 1.  The
pathwise identity
        |mu| tau = u + R + M,   u = log n,
        R = (tau - a) log 2 - a log(3/2) - log n   (plus one budget, >= 0)
        M = (log 3 / 2) (2a - tau)                   (centred parity sum)
holds for every n.  In the Bernoulli model R is the undershoot below 0 and
Wald's identities give E M = 0, E M^2 = sigma^2 E tau.  For each dyadic block
[2^k, 2^{k+1}) we compute the cumulants of tau and the moments of (u, R, M),
and we split M = M_head + M_tail, where M_head is the sum over the first k
steps (exactly fair coin flips by Terras) and M_tail the sum over the rest.
Output: data/ch3_decomp.npz  (one row per block)"""
import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
os.makedirs(BOOK + '/scripts/data', exist_ok=True)  # intermediate results that later scripts load
import numpy as np, numba as nb, sys, time
LOGN = int(sys.argv[1]) if len(sys.argv) > 1 else 30
N = 1 << LOGN
L2, L32, L3 = np.log(2.0), np.log(1.5), np.log(3.0)

@nb.njit(cache=True)
def census(N, tau, odd):
    for n in range(2, N):
        x = n; s = 0; a = 0
        while x >= n:
            if x & 1:
                x = (3*x + 1) >> 1; a += 1
            else:
                x >>= 1
            s += 1
        tau[n] = s + tau[x]; odd[n] = a + odd[x]

@nb.njit(cache=True)
def block_stats(tau, odd, k, L2, L32, L3):
    lo = 1 << k; hi = 1 << (k + 1); m = hi - lo
    # first pass: means
    st = 0.0; su = 0.0; sR = 0.0; sM = 0.0; sMh = 0.0; sMt = 0.0
    for n in range(lo, hi):
        t = tau[n]; a = odd[n]; u = np.log(n)
        R = (t - a)*L2 - a*L32 - u
        M = 0.5*L3*(2*a - t)
        # head: first k steps
        x = n; ah = 0
        for i in range(k):
            if x & 1:
                x = (3*x + 1) >> 1; ah += 1
            else:
                x >>= 1
        Mh = 0.5*L3*(2*ah - k)
        st += t; su += u; sR += R; sM += M; sMh += Mh; sMt += M - Mh
    et = st/m; eu = su/m; eR = sR/m; eM = sM/m; eMh = sMh/m; eMt = sMt/m
    c = np.zeros((6, 6)); k3 = 0.0; k4 = 0.0
    v = np.zeros(6)
    for n in range(lo, hi):
        t = tau[n]; a = odd[n]; u = np.log(n)
        R = (t - a)*L2 - a*L32 - u
        M = 0.5*L3*(2*a - t)
        x = n; ah = 0
        for i in range(k):
            if x & 1:
                x = (3*x + 1) >> 1; ah += 1
            else:
                x >>= 1
        Mh = 0.5*L3*(2*ah - k)
        v[0] = t - et; v[1] = u - eu; v[2] = R - eR; v[3] = M - eM
        v[4] = Mh - eMh; v[5] = (M - Mh) - eMt
        for i in range(6):
            for j in range(6):
                c[i, j] += v[i]*v[j]
        k3 += v[0]**3; k4 += v[0]**4
    c /= m; k3 /= m; k4 = k4/m - 3*c[0, 0]**2
    return et, eu, eR, eM, eMh, eMt, c, k3, k4

if __name__ == '__main__':
    t0 = time.time()
    tau = np.zeros(N, dtype=np.int16); odd = np.zeros(N, dtype=np.int16)
    census(N, tau, odd)
    print('census to 2^%d in %.0fs, max tau %d' % (LOGN, time.time() - t0, tau.max()), flush=True)
    np.save(BOOK + '/scripts/data/ch3_tau_small.npy', tau[:1 << 26])
    rows = []
    for k in range(3, LOGN):
        et, eu, eR, eM, eMh, eMt, c, k3, k4 = block_stats(tau, odd, k, L2, L32, L3)
        row = [k, et, eu, eR, eM, eMh, eMt, k3, k4] + list(c[np.triu_indices(6)])
        rows.append(row)
        print('k=%2d Etau %.3f var %.2f k3 %.0f | ER %.4f EM %.4f EMh %.4f EMt %.4f | VarR %.5f VarM %.3f VarMh %.3f VarMt %.3f CovMhMt %.4f CovRM %.4f  (%.0fs)'
              % (k, et, c[0, 0], k3, eR, eM, eMh, eMt, c[2, 2], c[3, 3], c[4, 4], c[5, 5], c[4, 5], c[2, 3], time.time() - t0), flush=True)
    np.savez(BOOK + '/scripts/data/ch3_decomp.npz', rows=np.array(rows))
