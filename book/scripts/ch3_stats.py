"""Chapter 3 computations: stopping time slope and dispersion, parity ACF,
excursion tail versus the Ville bound 1/K, standard map stopping slope."""
import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
import numpy as np
LOG2, LOG3 = np.log(2.0), np.log(3.0)
ds = 0.5*np.log(0.75)
sig2 = 0.5*LOG2**2 + 0.5*np.log(1.5)**2 - ds**2

def T(n): return n >> 1 if n % 2 == 0 else (3*n+1) >> 1
def T0(n): return n >> 1 if n % 2 == 0 else 3*n+1

def sigma_inf(n, f=T):
    s = 0
    while n != 1:
        n = f(n); s += 1
    return s

# A. stopping time slope (accelerated map), n in [2,10000]
ns = np.arange(2, 10001)
st = np.array([sigma_inf(int(n)) for n in ns])
L = np.log(ns)
g = (L @ st)/(L @ L)
print("A slope through origin, n<=1e4:", round(g, 4))
ns5 = np.arange(2, 5001); g5 = (np.log(ns5) @ st[:4999])/(np.log(ns5) @ np.log(ns5))
print("A slope n<=5000:", round(g5, 4))
print("model sigma2 %.6f  C %.4f  gamma0 %.5f  sig2/|ds|^3 %.3f" % (sig2, sig2/ds**2, 1/abs(ds), sig2/abs(ds)**3))
# dispersion by bins: larger range, random sample of n in [1e3, 1e7] log uniform
rng = np.random.default_rng(1)
for lo, hi in [(1e3, 2e3), (1e4, 2e4), (1e5, 2e5), (1e6, 2e6), (1e7, 2e7)]:
    sample = rng.integers(int(lo), int(hi), 4000)
    s = np.array([sigma_inf(int(n)) for n in sample]); l = np.log(sample)
    resid = s - l/abs(ds)
    mv = sig2/abs(ds)**3 * l.mean()
    print("A bin [%g,%g): mean st %.2f model mean %.2f | var %.1f model var %.1f ratio %.3f" %
          (lo, hi, s.mean(), (l/abs(ds)).mean(), s.var(), mv, s.var()/mv))
# standard map slope
st0 = np.array([sigma_inf(int(n), T0) for n in ns])
print("A standard map slope n<=1e4:", round((L @ st0)/(L @ L), 4), " model", round(3/np.log(4/3), 4))

# B. parity ACF, standard map
def acf(b, H):
    p = b.mean(); v = p*(1-p)
    return p, [(np.mean(b[:-h]*b[h:]) - p*p)/v for h in range(1, H+1)]
# B1 replicate small number scan (seed 0, start 27, restart odd < 10000), 4e5 steps
rng0 = np.random.default_rng(0)
bits, x = [], 27
while len(bits) < 400000:
    bits.append(x % 2); x = T0(x)
    if x == 1: x = 2*int(rng0.integers(1, 5000)) + 1
p, r = acf(np.array(bits, float), 8)
print("B1 small scan P(odd) %.4f rho %s" % (p, np.round(r, 4)))
# halving exponents K = gap - 1 between successive odd steps of the scan (lag 3 is (freq(K=2) - p)/(1 - p))
Kgap = np.diff(np.flatnonzero(np.array(bits) == 1)) - 1
print("B1 halving exponents: freq K=1 %.4f K=2 %.4f K=3 %.4f K>=4 %.4f" %
      ((Kgap == 1).mean(), (Kgap == 2).mean(), (Kgap == 3).mean(), (Kgap >= 4).mean()))
# B2 large random starts: n < 2^208 with bit 200 set (n mod 2^200 uniform), T0 steps 50 to 299 recorded, 4000 orbits
rngb = np.random.default_rng(2)
H = 10; num = np.zeros(H); cnt = 0; ones = 0; tot = 0; M = []
for i in range(4000):
    n = int.from_bytes(rngb.bytes(26), 'big') | (1 << 200)
    b = []
    for t in range(300):
        b.append(n & 1); n = T0(n)
    M.append(b[50:])   # drop a burn in of 50 steps
B = np.array(M, float)
p = B.mean(); v = p*(1-p)
rho = [ (np.mean(B[:, :-h]*B[:, h:]) - p*p)/v for h in range(1, H+1)]
print("B2 large starts P(odd) %.4f" % p)
print("B2 rho", np.round(rho, 4)); print("theory", np.round([(-0.5)**h for h in range(1, H+1)], 4))
np.save(BOOK + '/scripts/ch3_acf.npy', np.array([r + [0, 0], rho]) if False else np.array(rho))

# C. excursion tail: fraction of n in [2,1e6] with max_k T^k(n) >= K n
N = 10**6
tmax = np.zeros(N+1)
ratio = []
for n in range(2, N+1):
    m = n; x = n
    while x != 1:
        x = T(x)
        if x > m: m = x
    ratio.append(m/n)
ratio = np.array(ratio)
for K in [2, 4, 8, 16, 32, 64, 128, 256, 1024, 4096]:
    print("C K=%5d frac %.5f  bound %.5f  K*frac %.3f" % (K, (ratio >= K).mean(), 1/K, K*(ratio >= K).mean()))
lr = np.log(ratio*np.arange(2, N+1))/np.log(np.arange(2, N+1))
i = lr[1000:].argmax() + 1000
print("C max log t(n)/log n over 1000<n<=1e6: %.4f at n=%d" % (lr[i], i+2))
np.save(BOOK + '/scripts/ch3_ratio.npy', ratio)
