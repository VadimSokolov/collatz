import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
os.makedirs(BOOK + '/scripts/data', exist_ok=True)  # intermediate results that later scripts load
# Perfect split sampling for the energy tau (coefficient stopping time) under Haar measure.
# Level chain with the ideal weight Omega = 1/Z, truncated at depth K: m | x ~ omega on {0, ..., min(tau(x), K+1) - 1},
# x | m ~ exact alive word to depth m extended with fair bits until the word dies or reaches depth K (a residue
# alive at depth K is recorded as tau = K + 1, a censored value).  Check Polson and Scott: P(tau_{n+1} <= tau_n) = 1/2.
import numpy as np, math, json
K = 400; rng = np.random.default_rng(0)
N = [dict() for _ in range(K + 1)]; N[0][0] = 1
for k in range(1, K + 1):
    for m, c in N[k - 1].items():
        for mm in (m, m + 1):
            if 3 ** mm > 2 ** k: N[k][mm] = N[k].get(mm, 0) + c
Z = np.array([1.0] + [sum(N[k].values()) / 2 ** k for k in range(1, K + 1)])      # Z[k] = P(tau > k)
# suffix counts for exact conditional sampling to depth k, for every k (computed lazily)
cache = {}
def suffix(Kd):
    if Kd in cache: return cache[Kd]
    S = [dict() for _ in range(Kd + 1)]
    for m in N[Kd]: S[Kd][m] = 1
    for k in range(Kd - 1, -1, -1):
        for m in N[k]:
            S[k][m] = sum(S[k + 1].get(mm, 0) for mm in (m, m + 1) if 3 ** mm > 2 ** (k + 1))
    cache[Kd] = S; return S
def draw_tau_given_level(m):
    # word alive to depth m (exact), then fair bits until death; return tau
    S = suffix(m); ones = 0
    for k in range(m):
        c0 = S[k + 1].get(ones, 0) if 3 ** ones > 2 ** (k + 1) else 0
        c1 = S[k + 1].get(ones + 1, 0) if 3 ** (ones + 1) > 2 ** (k + 1) else 0
        ones += 1 if rng.random() < c1 / (c0 + c1) else 0
    k = m
    while k < K:
        k += 1; ones += rng.random() < 0.5
        if 3 ** ones < 2 ** k: return k
    return K + 1
# discrete ideal weight: omega_k proportional to 1/Z_k - 1/Z_{k-1} (k >= 1), omega_0 = 1
omega = np.concatenate([[1.0], 1 / Z[1:] - 1 / Z[:-1]])
def draw_level_given_tau(t):
    w = omega[:t]; w = w / w.sum(); return int(rng.choice(len(w), p=w))
t = draw_tau_given_level(0); ups = 0; same = 0; n = 20000; levels = []; taus = []
for i in range(n):
    m = draw_level_given_tau(t); t2 = draw_tau_given_level(m)
    ups += t2 > t; same += t2 == t; t = t2; levels.append(m); taus.append(t)
print("P(tau_{n+1} > tau_n) = %.4f  P(equal) = %.4f  P(down) = %.4f" % (ups / n, same / n, 1 - (ups + same) / n))
levels = np.array(levels); taus = np.array(taus)
print("level marginal: mean %.1f, max %d; tau visited up to %d; naive Monte Carlo would need %.1e draws to see tau > %d" % (
    levels.mean(), levels.max(), taus.max(), 1 / Z[min(taus.max(), K)], min(taus.max(), K)))
dl = np.diff(levels)
print("level moves between sweeps: up %.4f, down %.4f, equal %.4f" % ((dl > 0).mean(), (dl < 0).mean(), (dl == 0).mean()))
print("sweeps at the top level %d: %.4f (one in %.0f); sweeps returning tau > %d (recorded as %d): %.4f (one in %.1f)" % (
    K, (levels == K).mean(), 1 / max((levels == K).mean(), 1e-12), K, K + 1, (taus == K + 1).mean(), 1 / max((taus == K + 1).mean(), 1e-12)))
# exact stationary law of the truncated chain: level k has weight 1 (k = 0) or the hazard 1 - Z_k/Z_{k-1}
haz = 1 - Z[1:] / Z[:-1]; mass = 1 + haz.sum()
print("stationary: mass %.3f, mean level %.1f, P(level = %d) = %.4f (one in %.0f), P(tau > %d) = %.4f (one in %.1f)" % (
    mass, (np.arange(1, K + 1) * haz).sum() / mass, K, haz[-1] / mass, mass / haz[-1], K, 1 / mass, mass))
json.dump({'levels': levels.tolist(), 'taus': taus.tolist(), 'Z': Z.tolist()}, open(BOOK + '/scripts/data/ch10c_perfect.json', 'w'))
