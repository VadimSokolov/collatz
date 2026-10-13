import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
os.makedirs(BOOK + '/scripts/data', exist_ok=True)  # intermediate results that later scripts load
# Energy, level sets and nested sampling for the Collatz map.
# (1) Exact partition function of the energy E = tau (coefficient stopping time) under Haar prior:
#     Z(beta) = E[e^{beta tau}] = sum_k e^{beta k} P(tau = k), with P(tau > k) = 1 - F(k) from the DP.
#     Radius of convergence = the Chernoff exponent D.  Nested sampling evaluates exactly this kind
#     of integral from the level set volumes X_k = P(tau > k), here known exactly.
# (2) Nested sampling as a search for extremal integers: prior uniform on [2, 2^26), likelihood
#     L(n) = total stopping time of T (from the census table), MCMC moves = bit flips (2 adic
#     neighbours).  Compare the best found with the exact maximum.
import numpy as np, math, json, time
D = math.log(2) + (math.log(2)/math.log(3))*math.log(math.log(2)/math.log(3)) + (1 - math.log(2)/math.log(3))*math.log(1 - math.log(2)/math.log(3))  # log 2 - H(1/log2 3)
# ---- (1)
K = 400; N = {0: 1}; alive = []
for k in range(1, K + 1):
    new = {}
    for m, c in N.items():
        for mm in (m, m + 1):
            if 3 ** mm > 2 ** k: new[mm] = new.get(mm, 0) + c
    N = new; alive.append(sum(N.values()))
X = np.array([a / 2 ** k for k, a in zip(range(1, K + 1), alive)])      # X_k = P(tau > k)
P = np.diff(np.concatenate([[1.0], X])) * -1                              # P(tau = k), k = 1..K
print("P(tau=k) sums to", P.sum(), " tail X_400 =", X[-1])
for beta in (0.0, 0.01, 0.02, 0.03, 0.034, 0.0346, 0.035, 0.04):
    terms = np.exp(beta * np.arange(1, K + 1)) * P
    print(f"beta={beta:.4f}: partial Z over k<=400 = {terms.sum():.6f}; last term/first term ratio growth {terms[-1]/terms[-2]:.4f}")
print("exact decay of X_k: log X_400/400 =", math.log(X[-1]) / 400, " D =", D)
# ---- (2)
tau = np.load(BOOK + '/scripts/data/ch3_tau_small.npy')        # tau(n) of T for n < 2^26
M = len(tau); B = int(math.log2(M))
exact_max = int(tau[2:].max()); argmax = int(np.argmax(tau[2:]) + 2)
print("exact max total stopping time below 2^%d: %d at n = %d" % (B, exact_max, argmax))
rng = np.random.default_rng(3)
def nested(nlive=200, steps=20, iters=4000):
    live = rng.integers(2, M, nlive); L = tau[live].astype(float)
    best = (0, 0); hist = []
    for it in range(iters):
        i = int(np.argmin(L)); Lmin = L[i]; hist.append(Lmin)
        # replace with a prior sample constrained to L > Lmin: MCMC from a random other live point
        j = int(rng.integers(nlive));
        while j == i: j = int(rng.integers(nlive))
        x = int(live[j])
        for _ in range(steps):
            y = x ^ (1 << int(rng.integers(0, B)))        # flip one bit: a 2 adic neighbour
            if y >= 2 and tau[y] > Lmin: x = y
        live[i] = x; L[i] = tau[x]
        if L[i] > best[0]: best = (int(L[i]), x)
        if best[0] >= exact_max: return it, best, hist
    return iters, best, hist
t = time.time()
for seed in range(3):
    rng = np.random.default_rng(seed)
    it, best, hist = nested()
    print(f"nested sampling run {seed}: reached tau={best[0]} at n={best[1]} after {it} iterations ({(it+1)*20} likelihood calls; exact max {exact_max})")
print("time", time.time() - t)
# evidence estimate for L = e^{beta tau} by nested sampling vs exact sum over the table (beta=0.02)
beta = 0.02
rng = np.random.default_rng(1); it, best, hist = nested(nlive=400, steps=20, iters=6000)
hist = np.array(hist); nl = 400
Xk = np.exp(-np.arange(1, len(hist) + 1) / nl)         # shrinkage X_i = e^{-i/N}
w = -np.diff(np.concatenate([[1.0], Xk]))
Z_ns = np.sum(w * np.exp(beta * hist))
Z_exact = np.mean(np.exp(beta * tau[2:].astype(float)))
print(f"evidence E[e^{{beta tau}}] on [2,2^{B}): nested sampling {Z_ns:.4f}  exact {Z_exact:.4f}")
json.dump({'X': X.tolist(), 'P': P.tolist()}, open(BOOK + '/scripts/data/ch10b_levels.json', 'w'))
