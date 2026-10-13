"""Chapter 4: exact multi stride ascent densities and the passage time bounds in the stride model.

d_J = density among odd x of {S^J(x) > x} = sum_{K=J}^{floor(J lt)} C(K-1,J-1) 2^{-K}.
Exceptional finite sets E_J = {odd x : S^J(x) > x and K_J(x) > J lt}, contained in [1, X_J).
"""
import math
from fractions import Fraction
import numpy as np

LT = math.log2(3)
th = 1 / LT
D = th * math.log(2 * th) + (1 - th) * math.log(2 * (1 - th))

def syr(x):
    m = 3 * x + 1; k = 0
    while m % 2 == 0: m //= 2; k += 1
    return m, k

def v2(n):
    k = 0
    while n % 2 == 0: n //= 2; k += 1
    return k

# 1. stride exponent vectors define single classes mod 2^{K+1} (brute force K <= 14)
from collections import defaultdict
for J in (1, 2, 3, 4):
    M = 2 ** 16
    cls = defaultdict(set)
    for x in range(1, 2 ** 18, 2):
        y = x; ks = []
        for _ in range(J):
            y, k = syr(y); ks.append(k)
        K = sum(ks)
        if K + 1 <= 16:
            cls[tuple(ks)].add(x % 2 ** (K + 1))
    assert all(len(s) == 1 for s in cls.values())
print("stride vectors give single classes mod 2^(K+1) (checked J<=4, K<=15)")

# 2. exact densities and finite exceptional sets
def dJ(J):
    return sum(Fraction(math.comb(K - 1, J - 1), 2 ** K) for K in range(J, math.floor(J * LT) + 1))

print("J  d_J  float  X_J  E_J")
rows = []
for J in range(1, 21):
    Kp = math.floor(J * LT) + 1
    XJ = (1.5 ** J - 1) / (1 - 3 ** J / 2 ** Kp)
    E = []
    for x in range(1, int(XJ) + 2, 2):
        y = x; K = 0
        for _ in range(J):
            y, k = syr(y); K += k
        if y > x and K > J * LT:
            E.append(x)
    d = dJ(J)
    rows.append((J, d, XJ, E))
    print(J, d, "%.6f" % float(d), "%.1f" % XJ, E)

# 3. empirical check over odd x < 2^22
N = 2 ** 22
xs = np.arange(1, N, 2, dtype=np.int64)
y = xs.copy()
for J in range(1, 13):
    t = 3 * y + 1
    while True:
        ev = (t % 2 == 0)
        if not ev.any(): break
        t[ev] //= 2
    y = t
    print("J", J, "empirical frac S^J>x (odd x<2^22)", (y > xs).mean(), "exact", float(dJ(J)))

# 4. rate: d_J vs exp(-D lt J)
for J in (5, 10, 20, 40, 80, 160):
    d = float(dJ(J))
    print("J", J, "d_J %.4e" % d, "exp(-D lt J) %.4e" % math.exp(-D * LT * J), "-(1/J)log d_J %.5f" % (-math.log(d) / J))
print("D lt =", D * LT)
# Chernoff minimiser check
zs = np.linspace(0.01, 0.99, 100000)
phi = zs ** (1 - LT) / (2 - zs)
print("min phi", phi.min(), "at z", zs[phi.argmin()], "exp(-D lt)", math.exp(-D * LT), "2(c-1)/c", 2 * (LT - 1) / LT)

# 5. initial ascent run distribution among odd x < 2^20
cnt = defaultdict(int); tot = 0
for x in range(1, 2 ** 20, 2):
    cnt[v2(x + 1) - 1] += 1; tot += 1
print("initial run law", [(r, cnt[r] / tot, 2.0 ** (-(r + 1))) for r in range(6)])

# 6. strict ratio S^J(x)/x > (3/2)^J on x = -1 mod 2^{J+1}
for J in range(1, 9):
    for q in range(1, 400):
        x = 2 ** (J + 1) * q - 1
        y = x
        for _ in range(J): y, _k = syr(y)
        assert Fraction(y, x) > Fraction(3, 2) ** J
print("S^J(x)/x > (3/2)^J verified on the classes, J<=8")

# 7. sharp single stride descent bound: max over x = 1 mod 4, x >= 5, of S(x)/x is 7/9 at x = 9
best = max(((syr(x)[0] / x), x) for x in range(5, 10 ** 6, 4))
print("max S(x)/x over x=1 mod 4:", best, 7 / 9)

# 8. stride model passage time (prop:ch4:passage): by Wald's identity
#    E tau = (log x0 + E O)/mu with overshoot O = -log x_tau, and log 2 <= E O < 2 log 2,
#    so (log x0 + log 2)/mu <= E tau < (log x0 + 2 log 2)/mu; log x0/mu is only the leading term
mu = -math.log(0.75)
print("1/mu", 1 / mu, "log2/mu", math.log(2) / mu, "2 log2/mu", 2 * math.log(2) / mu)

def passage_exact(u0, tol=1e-15):
    """Exact E tau and E O by recursion on the law of the number K of halvings after n strides.

    f[K] = P(K halvings after n strides and tau > n); the walk is at u0 + n log 3 - K log 2,
    and the law of K after one more stride is g(K) = sum_{k>=1} f(K-k) 2^{-k} = (g(K-1) + f(K-1))/2.
    """
    f = {0: 1.0}; n = 0; Et = 0.0; EO = 0.0
    while True:
        alive = sum(f.values())
        Et += alive                      # E tau = sum_{n>=0} P(tau > n)
        if alive < tol: break
        g = {}; prev = 0.0
        for K in range(1, max(f) + 81):  # mass beyond max(f)+80 is below 2^{-80}
            prev = (prev + f.get(K - 1, 0.0)) / 2
            u = u0 + (n + 1) * math.log(3) - K * math.log(2)
            if u > 0:
                g[K] = prev
            else:
                EO += prev * (-u)        # absorbed at stride n+1 with overshoot -u
        f = g; n += 1
    return Et, EO

for u0 in (2.0, 5.0, 10.0, 20.0):
    Et, EO = passage_exact(u0)
    print("exact: log x0 %.1f  E tau %.4f  E O %.4f  (log x0 + E O)/mu %.4f" % (u0, Et, EO, (u0 + EO) / mu))

# simulation check of the exact values
rng = np.random.default_rng(7)
for u0 in (2.0, 5.0, 10.0, 20.0):
    n = 200000
    u = np.full(n, u0); tau = np.zeros(n); alive = np.ones(n, bool)
    over = np.zeros(n)
    while alive.any():
        k = rng.geometric(0.5, size=alive.sum())
        u[alive] += math.log(3) - k * math.log(2)
        tau[alive] += 1
        newly = alive.copy(); newly[alive] = u[alive] <= 0
        over[newly] = -u[newly]
        alive &= ~newly
    print("simulation: log x0 %.1f  E tau %.4f (SE %.4f)  leading %.4f  lower %.4f  upper %.4f  E O %.4f  Wald %.4f" %
          (u0, tau.mean(), tau.std() / math.sqrt(n), u0 / mu, (u0 + math.log(2)) / mu,
           (u0 + 2 * math.log(2)) / mu, over.mean(), (u0 + over.mean()) / mu))
