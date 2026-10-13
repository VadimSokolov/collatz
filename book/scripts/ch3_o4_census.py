"""Chapter 3, second order structure of stopping times (problem O4).
Census of the total stopping time sigma_inf(n) of the accelerated map T for
every n < 2^LOGN, by memoisation in numba, and blockwise moments over the
dyadic blocks [2^k, 2^{k+1}).  Also the first passage time below n (the
stopping time sigma(n)), and the landing value at the first passage below n.
Output: data/ch3_census.npz"""
import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
os.makedirs(BOOK + '/scripts/data', exist_ok=True)  # intermediate results that later scripts load
import numpy as np, numba as nb, sys, time
LOGN = int(sys.argv[1]) if len(sys.argv) > 1 else 30
N = 1 << LOGN
LIM = np.uint64((1 << 63) // 3)

@nb.njit(cache=True)
def census(N, memo):
    over = 0
    for n in range(2, N):
        x = np.uint64(n); s = 0
        while x >= n:
            if x & np.uint64(1):
                if x > LIM:
                    over += 1; break
                x = (np.uint64(3)*x + np.uint64(1)) >> np.uint64(1)
            else:
                x = x >> np.uint64(1)
            s += 1
        memo[n] = s + memo[x]
    return over

@nb.njit(cache=True)
def block_moments(memo, k):
    lo = 1 << k; hi = 1 << (k+1); m = hi - lo
    s1 = 0.0; l1 = 0.0
    for n in range(lo, hi):
        s1 += memo[n]; l1 += np.log(n)
    mean = s1/m; lmean = l1/m
    c2 = 0.0; c3 = 0.0; c4 = 0.0; cl = 0.0; vl = 0.0
    for n in range(lo, hi):
        d = memo[n]-mean; c2 += d*d; c3 += d*d*d; c4 += d*d*d*d
        e = np.log(n) - lmean; cl += d*e; vl += e*e
    return mean, c2/m, c3/m, c4/m, lmean, cl/m, vl/m

@nb.njit(cache=True)
def stop_moments(k):
    # sigma(n): first j with T^j(n) < n, for n in block k (n>=2)
    lo = max(1 << k, 2); hi = 1 << (k+1)
    s1 = 0.0; s2 = 0.0; m = 0
    for n in range(lo, hi):
        x = np.uint64(n); s = 0
        while x >= n:
            if x & np.uint64(1):
                x = (np.uint64(3)*x + np.uint64(1)) >> np.uint64(1)
            else:
                x = x >> np.uint64(1)
            s += 1
        s1 += s; s2 += s*s; m += 1
    return s1/m, s2/m - (s1/m)**2

t0 = time.time()
memo = np.zeros(N, dtype=np.int16)
over = census(N, memo)
print('census to 2^%d done in %.0fs, overflow %d, max %d' % (LOGN, time.time()-t0, over, memo.max()))
rows = []
for k in range(3, LOGN):
    mean, v, c3, c4, lm, cl, vl = block_moments(memo, k)
    sm, sv = stop_moments(k) if k <= 26 else (np.nan, np.nan)
    rows.append((k, mean, v, c3/v**1.5, c4/v**2, lm, cl, vl, sm, sv))
    print('k=%2d mean %8.3f var %9.3f skew %.4f kurt %.4f E log n %.4f stop mean %.3f var %.3f'
          % (k, mean, v, c3/v**1.5, c4/v**2, lm, sm, sv))
np.savez(BOOK + '/scripts/data/ch3_census.npz', rows=np.array(rows))
np.save(BOOK + '/scripts/data/ch3_memo_small.npy', memo[:1 << 24])
