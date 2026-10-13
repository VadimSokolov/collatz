"""Chapter 2: exact distribution of coefficient stopping times over residues mod 2^k.

N_k(m) = number of words w in {0,1}^k with 3^{m_j(w)} > 2^j for all 1<=j<=k and m_k(w)=m.
Recursion N_k(m) = (N_{k-1}(m) + N_{k-1}(m-1)) * [3^m > 2^k], N_0(0)=1.
R_k = number of residues mod 2^k with tau = k = N_{k-1}(m*) where 2^{k-1} <= 3^{m*} < 2^k
(equality on the left only for k = 1, m* = 0).
Checks the ballot sandwich P(X >= theta k)/k <= 1-F(k) <= P(X >= theta k), X ~ Bin(k,1/2),
and the explicit bound 1-F(k) >= exp(-Dk-1)/(k(k+1)).
"""
import math
from fractions import Fraction
from collections import defaultdict

KMAX = 400
th = math.log(2) / math.log(3)
D = th * math.log(2 * th) + (1 - th) * math.log(2 * (1 - th))
print("theta", th, "D", D)

def T(n): return n // 2 if n % 2 == 0 else (3 * n + 1) // 2

N = [{0: 1}]
R = [None]
for k in range(1, KMAX + 1):
    prev = N[-1]
    new = defaultdict(int)
    for m, c in prev.items():
        for b in (0, 1):
            mm = m + b
            if 3 ** mm > 2 ** k:
                new[mm] += c
    N.append(dict(new))
    # residues with tau = k: survivors at k-1 with 2^{k-1} <= 3^m < 2^k, followed by a 0
    r = sum(c for m, c in prev.items() if 3 ** m < 2 ** k)
    ms = [m for m in prev if 3 ** m < 2 ** k]
    R.append(r)
    assert len(ms) <= 1

surv = [1] + [sum(N[k].values()) for k in range(1, KMAX + 1)]
# consistency: surv_k = 2 surv_{k-1} - R_k
for k in range(1, KMAX + 1):
    assert surv[k] == 2 * surv[k - 1] - R[k], k
# R_k > 0 iff k = floor(m log2 3) + 1
allowed = set(math.floor(m * math.log2(3)) + 1 for m in range(0, KMAX))
for k in range(1, KMAX + 1):
    assert (R[k] > 0) == (k in allowed), k
print("R_k > 0 exactly at k = floor(m lt)+1, k <=", KMAX)

# brute force for k <= 20
def tau_res(r, k):
    m = 0; x = r
    for j in range(1, k + 1):
        m += x % 2; x = T(x)
        if 3 ** m < 2 ** j: return j
    return None
for k in range(1, 21):
    cnt = sum(1 for r in range(2 ** k) if tau_res(r, k) == k)
    assert cnt == R[k], (k, cnt, R[k])
print("brute force R_k agrees for k <= 20")

# table to k = 40
print("k R_k survivors 1-F(k)")
for k in range(1, 41):
    print(k, R[k], surv[k], "%.6f" % (surv[k] / 2 ** k))

# sandwich
def binom_tail(k):
    m0 = math.ceil(th * k)
    return Fraction(sum(math.comb(k, i) for i in range(m0, k + 1)), 2 ** k)
worst_low = 1e9; worst_exp = 1e9
for k in range(1, KMAX + 1):
    q = Fraction(surv[k], 2 ** k)
    bt = binom_tail(k)
    assert bt / k <= q <= bt, k
    lowb = math.exp(-D * k - 1) / (k * (k + 1))
    assert float(q) >= lowb and float(q) <= math.exp(-D * k), k
    worst_low = min(worst_low, float(q * k / bt))
print("sandwich and explicit bounds hold for k <=", KMAX, "; min of (1-F)k/tail =", worst_low)
for k in (20, 40, 100, 200, 400):
    q = surv[k] / 2 ** k; bt = float(binom_tail(k))
    print(k, "1-F", "%.4e" % q, "tail", "%.4e" % bt, "chernoff", "%.4e" % math.exp(-D * k),
          "-(1/k)log(1-F)", "%.5f" % (-math.log(q) / k),
          "k^1.5 e^{Dk}(1-F)", "%.4f" % (k ** 1.5 * math.exp(D * k) * q),
          "(1-F)/tail", "%.4f" % (q / bt))
print("sum_k R_k/2^k up to 400 =", 1 - surv[KMAX] / 2 ** KMAX)

# sigma = tau for 2 <= n <= 10^6, and max sigma
Nmax = 10 ** 6
bad = 0; best = (0, 0)
for n in range(2, Nmax + 1):
    x = n; k = 0; m = 0; tau = None
    while True:
        m += x % 2; x = T(x); k += 1
        if tau is None and 3 ** m < 2 ** k: tau = k
        if x < n: break
    if tau != k: bad += 1
    if k > best[0]: best = (k, n)
print("sigma != tau count for 2..1e6:", bad, "max sigma", best)

# local exponent of the ratio (1-F(k))e^{Dk} between two values of k, by the recursion to k = 800
st = {0: 1}; surv800 = [1]
for k in range(1, 801):
    new = defaultdict(int)
    for m, c in st.items():
        for b in (0, 1):
            if 3 ** (m + b) > 2 ** k: new[m + b] += c
    st = new; surv800.append(sum(st.values()))
def log_ratio(k): return math.log(surv800[k]) - k * math.log(2) + D * k
for k1, k2 in ((20, 100), (100, 200), (200, 400), (400, 800)):
    print("local exponent of (1-F(k))e^{Dk}, k =", k1, "to", k2, ":",
          "%.3f" % ((log_ratio(k2) - log_ratio(k1)) / math.log(k2 / k1)))
