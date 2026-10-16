"""Out of sample test of the O4 fit: sample n uniform in [2^k, 2^{k+1}) for k = 32 .. 48, compute the total
stopping time of T exactly (iterate until below 2^26, then read the census table), and compare the
observed variance, third cumulant and skewness with the fit kappa_j = A_j E[log n] - D_j obtained from
the census blocks 20 <= k <= 29.  Orbits whose values would overflow 63 bits are finished in Python."""
import numpy as np, numba as nb, json, math, os
HERE = os.path.dirname(os.path.abspath(__file__))
tau = np.load(os.path.join(HERE, 'data', 'ch3_tau_small.npy')).astype(np.int64)
fit = json.load(open(os.path.join(HERE, 'data', 'ch16_o4.json')))['fit']
LIM = 1 << 26
@nb.njit(cache=True)
def stop_times(ns, tau, LIM):
    out = np.empty(ns.shape[0], np.int64); big = 1 << 61
    for i in range(ns.shape[0]):
        x = ns[i]; s = 0; ok = True
        while x >= LIM:
            if x & 1:
                if x > big: ok = False; break
                x = (3 * x + 1) >> 1
            else: x >>= 1
            s += 1
        out[i] = s + tau[x] if ok else -1
    return out
def slow(n):
    x = n; s = 0
    while x >= LIM:
        x = (3 * x + 1) >> 1 if x & 1 else x >> 1; s += 1
    return s + int(tau[x])
rng = np.random.default_rng(1601); res = []
for k in (32, 36, 40, 44, 48):
    ns = rng.integers(1 << k, 1 << (k + 1), 2_000_000, dtype=np.int64)
    t = stop_times(ns, tau, LIM); bad = np.nonzero(t < 0)[0]
    for i in bad: t[i] = slow(int(ns[i]))
    t = t.astype(float); u = np.log(ns.astype(float)).mean(); c = t - t.mean()
    v, k3 = c.var(), (c ** 3).mean()
    pv = fit['k2']['slope'] * u - fit['k2']['D']; pk3 = fit['k3']['slope'] * u - fit['k3']['D']
    # bootstrap standard error of the skewness
    sk = k3 / v ** 1.5; bs = []
    for b in range(20):
        idx = rng.integers(0, t.size, t.size); cb = t[idx] - t[idx].mean(); bs.append((cb ** 3).mean() / cb.var() ** 1.5)
    res.append(dict(k=k, mean=float(t.mean()), var=v, k3=k3, skew=sk, skew_se=float(np.std(bs)), pred_var=pv, pred_k3=pk3,
                    pred_skew=pk3 / pv ** 1.5, model_skew=(3 * (math.log(3) / 2) ** 4 / abs(math.log(math.sqrt(3) / 2)) ** 5 * u) /
                    ((math.log(3) / 2) ** 2 / abs(math.log(math.sqrt(3) / 2)) ** 3 * u) ** 1.5, overflow=int(bad.size)))
    r = res[-1]
    print(f"k={k}: var {v:.1f} (pred {pv:.1f}), k3 {k3:.0f} (pred {pk3:.0f}), skew {sk:.4f} +- {r['skew_se']:.4f} (pred {r['pred_skew']:.4f}, model {r['model_skew']:.4f}); overflow {bad.size}", flush=True)
json.dump(res, open(os.path.join(HERE, 'data', 'ch16_o4_predict.json'), 'w'), default=float)
