import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
os.makedirs(BOOK + '/scripts/data', exist_ok=True)  # intermediate results that later scripts load
# Level set Monte Carlo for the Collatz map: nested sampling (Skilling), split sampling (Birge,
# Chang and Polson) and vertical likelihood weighted slice sampling (Polson and Scott), applied to
#   prior   pi = uniform on [2, 2^B)  (Haar measure truncated),
#   L(x)    = total stopping time sigma(x) of the accelerated map T (census table, exact),
# estimating the rare event probabilities Z(m) = P(sigma > m), against the exact values.
# Also: perfect (exact) split sampling for the energy tau (coefficient stopping time), where the
# conditional pi(x | tau > m) can be sampled exactly by dynamic programming.
import numpy as np, math, json, time, sys
B = 26
tau = np.load(BOOK + '/scripts/data/ch3_tau_small.npy'); M = 1 << B
sig = tau[2:].astype(np.int64); cnt = np.bincount(sig, minlength=602); exactZ = (cnt.sum() - np.cumsum(cnt))[:601] / cnt.sum()
L2, L3 = math.log(2), math.log(3); theta = L2 / L3
print("exact: P(sigma>300) %.3e  P(sigma>400) %.3e  P(sigma>500) %.3e  max %d" % (exactZ[300], exactZ[400], exactZ[500], sig.max()))

# ---------- perfect split sampling for E = tau: exact conditional sampling of alive words ----------
def alive_counts(K):
    # N[k][m] = number of alive words of length k with m ones (3^j > 2^... prefix condition)
    N = [dict() for _ in range(K + 1)]; N[0][0] = 1
    for k in range(1, K + 1):
        for m, c in N[k - 1].items():
            for mm in (m, m + 1):
                if 3 ** mm > 2 ** k: N[k][mm] = N[k][mm] + c if mm in N[k] else c
    return N
def sample_alive_word(N, K, rng):
    # sample uniformly among alive words of length K by backward DP: suffix counts
    # S[k][m] = number of ways to extend a prefix of length k with m ones to an alive word of length K
    S = [dict() for _ in range(K + 1)]
    for m in N[K]: S[K][m] = 1
    for k in range(K - 1, -1, -1):
        for m in N[k]:
            S[k][m] = sum(S[k + 1].get(mm, 0) for mm in (m, m + 1) if 3 ** mm > 2 ** (k + 1))
    w = []; m = 0
    for k in range(K):
        c0 = S[k + 1].get(m, 0) if 3 ** m > 2 ** (k + 1) else 0
        c1 = S[k + 1].get(m + 1, 0) if 3 ** (m + 1) > 2 ** (k + 1) else 0
        b = 1 if rng.random() < c1 / (c0 + c1) else 0
        w.append(b); m += b
    return w
K = 60; N = alive_counts(K); rng = np.random.default_rng(0)
ws = [sample_alive_word(N, K, rng) for _ in range(3000)]
ones = np.array([sum(w) for w in ws]); tot = sum(N[K].values())
print("perfect split sampling at level %d: alive words %d (volume %.3e); sampled mean ones %.3f, exact %.3f" % (
    K, tot, tot / 2 ** K, ones.mean(), sum(m * c for m, c in N[K].items()) / tot))
# Terras: the word determines the residue class mod 2^K; recover a residue and check tau >= K
def residue_from_word(w):
    # find r mod 2^K whose parity vector is w (bijection): solve bit by bit
    r = 0
    for k in range(len(w)):
        x = r; ok = True
        for j in range(k + 1):
            if (x & 1) != w[j]: ok = False; break
            x = (3 * x + 1) // 2 if x & 1 else x // 2
        if not ok: r += 1 << k
    return r
r = residue_from_word(ws[0]); x = r; 
for j in range(K):
    assert (x & 1) == ws[0][j]; x = (3 * x + 1) // 2 if x & 1 else x // 2
print("word to residue bijection checked for a sample; residue", r)

# ---------- level set samplers for L = sigma on [2, 2^B) ----------
def mcmc_level(x, Lmin, rng, steps):
    for _ in range(steps):
        y = x
        for _ in range(int(rng.integers(1, 4))): y ^= 1 << int(rng.integers(0, B))
        if y >= 2 and tau[y] > Lmin: x = y
    return x
def nested(rng, nlive, iters, steps):
    live = rng.integers(2, M, nlive); L = tau[live].astype(float); Lm = []
    for it in range(iters):
        i = int(np.argmin(L)); Lmin = L[i]; Lm.append(Lmin)
        j = i
        while j == i: j = int(rng.integers(nlive))
        x = mcmc_level(int(live[j]), Lmin, rng, steps); live[i] = x; L[i] = tau[x]
    Lm = np.array(Lm); X = np.exp(-np.arange(1, iters + 1) / nlive)
    # Z(m) estimate: prior volume above m = X at the first iteration whose Lmin exceeds m
    def Zhat(m):
        k = np.searchsorted(Lm, m, side='right')      # Lm nondecreasing
        return X[k - 1] if k > 0 else 1.0
    Zhat.trace = dict(calls=iters * steps + nlive, Lfinal=float(Lm[-1]), live=live.copy())
    return Zhat, (iters * steps + nlive)
def split(rng, nsamp, steps, rho=math.exp(-1), mmax=560):
    # sequential levels at the (1-rho) quantile; product estimator Z_t = prod of conditional fractions
    levels = [0.0]; Zs = [1.0]; x = int(rng.integers(2, M)); calls = 0; last = None
    while levels[-1] < mmax:
        Lmin = levels[-1]; xs = []; pts = []
        while tau[x] <= Lmin: x = int(rng.integers(2, M)); calls += 1
        for _ in range(nsamp):
            x = mcmc_level(x, Lmin, rng, steps); xs.append(tau[x]); pts.append(x); calls += steps
        xs = np.array(xs, float); mt = np.quantile(xs, 1 - rho)
        if mt <= Lmin: mt = Lmin + 1
        frac = np.mean(xs > mt)
        if frac == 0: last = dict(level=float(Lmin), pts=sorted(set(pts)), top=float(xs.max()), attop=float(np.mean(xs == xs.max()))); break
        levels.append(mt); Zs.append(Zs[-1] * frac)
        while tau[x] <= mt: x = mcmc_level(x, Lmin, rng, 5); calls += 5
    levels = np.array(levels); Zs = np.array(Zs)
    trace = dict(calls=calls, stages=len(levels) - 1, top=float(levels[-1]), last=last)
    def Zhat(m):
        k = np.searchsorted(levels, m, side='right') - 1
        # interpolate piecewise exponentially between levels
        if k >= len(levels) - 1: return Zs[-1]
        f = (m - levels[k]) / (levels[k + 1] - levels[k])
        return Zs[k] * (Zs[k + 1] / Zs[k]) ** f
    Zhat.trace = trace
    return Zhat, calls
def vertical(rng, nsweep, steps, eta=1e-7, grid=None):
    # weighted slice sampling with W(u) = 1 / max(eta, Zhat(u)), Zhat learned on the fly from the
    # level marginal; importance weights 1/W(L(x)); Rao Blackwellised Z(m) estimate (split sampling eq.)
    if grid is None: grid = np.arange(0, 601)
    Zcur = np.exp(-0.02 * grid)                      # initial guess, refined between passes
    x = int(rng.integers(2, M)); calls = 0; samples = []
    for sweep in range(nsweep):
        Lx = tau[x]
        # draw u | x: with W = 1/max(eta,Z): T ~ Unif(0, 1/max(eta, Z(Lx))); u = 0 if T<=1 else Zinv(1/T)
        Tv = rng.uniform(0, 1 / max(eta, Zcur[min(Lx, 600)]))
        if Tv <= 1: u = 0.0
        else:
            target = 1 / Tv; k = np.searchsorted(-Zcur, -target); u = float(min(k, 600))
        x = mcmc_level(x, u - 1e-9, rng, steps); calls += steps; samples.append(tau[x])
        if sweep % 2000 == 1999:                       # refine Z from the samples so far
            s = np.array(samples); w = np.array([max(eta, Zcur[min(v, 600)]) for v in s])  # 1/W
            est = np.array([np.sum(w * (s > m)) / np.sum(w) for m in grid])
            est = np.maximum.accumulate(est[::-1])[::-1]; est[est <= 0] = est[est > 0].min()
            Zcur = est
    s = np.array(samples); w = np.array([max(eta, Zcur[min(v, 600)]) for v in s])
    def Zhat(m): return np.sum(w * (s > m)) / np.sum(w)
    Zhat.trace = dict(calls=calls, maxL=int(s.max()))
    return Zhat, calls
ms = [100, 200, 300, 400, 450, 500, 550]
res = {'exact': [exactZ[m] for m in ms]}
budget = 200000; traces = {}
for name, fn in (('naive', None), ('nested', nested), ('split', split), ('vertical', vertical)):
    ests = []
    for rep in range(6):
        rng = np.random.default_rng(100 + rep)
        if name == 'naive':
            xs = rng.integers(2, M, budget); Zh = lambda m, xs=xs: np.mean(tau[xs] > m)
            Zh.trace = dict(calls=budget, hits={m: int(np.sum(tau[xs] > m)) for m in (300, 400, 450, 500)})
        elif name == 'nested': Zh, calls = nested(rng, nlive=200, iters=1600, steps=120)
        elif name == 'split': Zh, calls = split(rng, nsamp=1200, steps=10)
        else: Zh, calls = vertical(rng, nsweep=12000, steps=16)
        ests.append([Zh(m) for m in ms]); traces.setdefault(name, []).append(Zh.trace)
    ests = np.array(ests, float); res[name] = ests.tolist()
    rel = np.sqrt(np.mean((np.log(np.maximum(ests, 1e-12)) - np.log(np.array(res['exact'])))**2, axis=0))
    print(name, 'RMSE of log Z at m =', ms, ':', np.round(rel, 3), ' mean estimates', ["%.2e" % v for v in ests.mean(0)])
print('exact', ["%.2e" % v for v in res['exact']])
print('medians', {k: ["%.3g" % v for v in np.median(np.array(v), axis=0)] for k, v in res.items() if k != 'exact'})
# ---------- trace: budgets, naive hits, and where the level set samplers stop ----------
def Tm(x): return (3 * x + 1) // 2 if x & 1 else x // 2
def key(n, r):                                     # orbit point of n < 2^B with remaining stopping time r
    x = int(n)
    for _ in range(int(tau[n]) - r): x = Tm(x)
    return x
for name in ('naive', 'nested', 'split', 'vertical'):
    print(name, 'energy evaluations per run:', [t['calls'] for t in traces[name]])
print('naive hits per run above m:', [t['hits'] for t in traces['naive']])
print('split: stages per run', [t['stages'] for t in traces['split']], ' top level', [t['top'] for t in traces['split']])
print('vertical: deepest energy visited per run', [t['maxL'] for t in traces['vertical']])
for t in traces['split']:                          # the stage at which split sampling stopped
    la = t['last']; lv = int(la['level']); pts = la['pts']
    print('split: stopped in the stage above level %d: %d distinct integers visited in %d classes (remaining time %d),'
          ' top energy %d held by %.2f of the samples' % (lv, len(pts), len({key(x, lv - 100) for x in pts}), lv - 100, la['top'], la['attop']))
for t in traces['nested']:
    mf = int(t['Lfinal']); r = mf - 100; live = t['live']
    S = np.nonzero(tau[2:] > mf)[0] + 2; ks = {}
    for n in S: ks.setdefault(key(n, r), []).append(int(tau[n]))
    occ = {key(x, r) for x in live if tau[x] > mf}
    print('nested: final threshold %d, deepest live point %d, %d live points above it in %d of the %d classes of {sigma > %d}'
          ' (classes at remaining time %d); deepest member of an occupied class %d' % (
          mf, int(tau[live].max()), int(np.sum(tau[live] > mf)), len(occ), len(ks), mf, r, max([max(ks[k]) for k in occ] + [mf])))
json.dump({'ms': ms, 'res': res, 'exactZ': exactZ.tolist()}, open(BOOK + '/scripts/data/ch10c_mc.json', 'w'))
