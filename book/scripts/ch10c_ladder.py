import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
os.makedirs(BOOK + '/scripts/data', exist_ok=True)  # intermediate results that later scripts load
# The exact staircase Z(k) = P(tau > k) of the coefficient stopping time under Haar measure, far
# beyond the depth 400 of the perfect sampler: the dynamic programme of the exact stopping counts,
# in floating point with the factor e^{-Dk} taken out (an exact integer check to depth 400 first).
# Prints: E[tau]; 1/Z(400); the sum of the hazards (the mass of the augmented law truncated at K);
# the mean hazard on the stopping levels k = floor(j log2 3) + 1; fits of the polynomial factor
# Z(k) e^{Dk} = C k^{-p} over several ranges; and the partition function E[e^{D tau}] at the
# boundary beta = D, with the Wald identity bound 2^{theta*}.
import numpy as np, math, json
L2, L3 = math.log(2), math.log(3); th = L2 / L3                   # th = 1/log2 3
D = L2 + th * math.log(th) + (1 - th) * math.log(1 - th)          # log 2 - H(1/log2 3), nats
print("D = %.6f" % D)
# exact integers to depth 400
N = {0: 1}; Zex = [1.0]
for k in range(1, 401):
    new = {}
    for j, c in N.items():
        for jj in (j, j + 1):
            if 3 ** jj > 2 ** k: new[jj] = new.get(jj, 0) + c
    N = new; Zex.append(sum(N.values()) / 2 ** k)
print("exact: Z(400) = %.6e, 1/Z(400) = %.4e" % (Zex[400], 1 / Zex[400]))
# floating point to depth K, q[j] = P(alive at step k with j ones) e^{Dk}
K = 20000; eD = math.exp(D)
q = np.zeros(K + 2); q[0] = 1.0; Zq = np.empty(K + 1); Zq[0] = 1.0
for k in range(1, K + 1):
    new = np.empty_like(q); new[0] = q[0]; new[1:] = q[1:] + q[:-1]; new *= 0.5 * eD
    new[:math.floor(k * th) + 1] = 0.0                            # alive: 3^j > 2^k, j > k th
    q = new; Zq[k] = q.sum()
k = np.arange(K + 1); logZ = np.log(Zq) - D * k                  # log Z(k)
print("float vs exact at k = 400: relative difference %.2e" % abs(math.exp(logZ[400]) / Zex[400] - 1))
Z = np.exp(logZ)
print("E[tau] = sum_k Z(k) = %.5f" % Z.sum())
h = np.concatenate([[0.0], 1 - np.exp(logZ[1:] - logZ[:-1])])    # hazard h_k, k >= 1
stop, p3 = [], 1                                                  # stopping levels k = bit length of 3^j
while p3.bit_length() <= K: stop.append(p3.bit_length()); p3 *= 3   # = floor(j log2 3) + 1, j >= 0
stop = np.array(stop); off = np.setdiff1d(np.arange(1, K + 1), stop)
print("%d stopping levels up to %d; largest |hazard| off them %.1e (rounding)" % (len(stop), K, np.abs(h[off]).max()))
for KK in (400, 1000, 6000, 20000):
    print("K = %5d: sum of hazards %.3f, mass of the truncated augmented law %.3f, stationary P(tau > K) = %.4f"
          % (KK, h[1:KK + 1].sum(), 1 + h[1:KK + 1].sum(), 1 / (1 + h[1:KK + 1].sum())))
for a, b in ((100, 400), (100, 6000), (1000, 6000), (10000, 20000)):
    s = stop[(stop > a) & (stop <= b)]
    print("mean hazard on the stopping levels %d < k <= %d: %.4f (%d levels)" % (a, b, h[s].mean(), len(s)))
print("pure exponential with one stopping level every log2 3 steps: 1 - e^{-D log2 3} = %.4f" % (1 - math.exp(-D * L3 / L2)))
y = logZ + D * k
for a, b in ((100, 400), (400, 1000), (1000, 2000), (1000, 6000), (6000, 20000)):
    for name, idx in (("all k", np.arange(a, b + 1)), ("stopping levels", stop[(stop >= a) & (stop <= b)])):
        A = np.vstack([np.ones(len(idx)), -np.log(idx)]).T; c = np.linalg.lstsq(A, y[idx], rcond=None)[0]
        print("fit Z(k) = C k^-p e^{-Dk} over %5d <= k <= %5d (%s): p = %.3f, C = %.3f" % (a, b, name, c[1], math.exp(c[0])))
terms = np.concatenate([[0.0], eD * Zq[:-1] - Zq[1:]])            # (Z(k-1) - Z(k)) e^{Dk}
S = np.cumsum(terms)
for KK in (100, 200, 400, 800, 1000, 2000, 6000, 10000, 20000):
    print("partition function at beta = D, partial sum to k = %5d: %.5f" % (KK, S[KK]))
c32 = (S[20000] - S[10000]) / (2 * (1 / math.sqrt(10000) - 1 / math.sqrt(20000)))
print("tail fit c k^{-3/2}: c = %.4f, extrapolated limit %.4f" % (c32, S[20000] + 2 * c32 / math.sqrt(20000)))
tstar = math.log(L2 / math.log(1.5)) / L3
rho = ((1.5) ** tstar + 0.5 ** tstar) / 2
print("theta* = %.4f, min E[xi^theta] = %.6f, -log = %.6f (= D), Wald bound 2^theta* = %.4f" % (tstar, rho, -math.log(rho), 2 ** tstar))
json.dump({'K': K, 'logZ': logZ.tolist()[:6001], 'D': D}, open(BOOK + '/scripts/data/ch10c_ladder.json', 'w'))
