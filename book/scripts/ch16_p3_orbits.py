"""P3: single orbits from very large starting values.  Syracuse orbits from random odd y_0 near 2^100000
(four orbits) and 2^300000 (one orbit).  Along each orbit:
  (1) residues y_j mod 3^a, a = 5..8: total variation distance of the empirical law from the stationary
      law pi_a of the induced chain, against the same distance for i.i.d. samples from pi_a of equal length;
  (2) consecutive residues mod 3^4: chi squared of the transition counts against the exact transition law
      of the induced chain (degrees of freedom from the nonzero cells);
  (3) halving exponents k_j: frequencies against Geometric(1/2), autocorrelations at lags 1..5, and the
      running mean against 2;
  (4) low binary digits y_j mod 2^12 against the uniform law on odd classes (chi squared).
Output data/ch16_p3.json."""
import numpy as np, random, math, json, os, time
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__))

def chain(a):
    M = 3 ** a; cls = np.array([r for r in range(M) if r % 3]); idx = -np.ones(M, int); idx[cls] = np.arange(len(cls))
    n = len(cls); per = 2 * 3 ** (a - 1); inv2 = pow(2, -1, M)
    P = np.zeros((n, n))
    pw = np.array([pow(inv2, K, M) for K in range(1, per + 1)], dtype=np.int64)
    wts = 2.0 ** -np.arange(1, per + 1) / (1 - 2.0 ** -per)
    for r in cls:
        v = ((3 * r + 1) % M * pw) % M
        np.add.at(P[:, idx[r]], idx[v], wts)
    pi = np.ones(n) / n
    for _ in range(4 * a + 10): pi = P @ pi
    return cls, idx, P, pi / pi.sum()

def orbit(y, a_max=8, b=12):
    M = 3 ** a_max; B = (1 << b) - 1; res = []; ks = []; low = []
    while y > 1:
        y = 3 * y + 1
        k = (y & -y).bit_length() - 1; y >>= k
        ks.append(k); res.append(y % M); low.append(y & B)
    return np.array(res, np.int64), np.array(ks, np.int64), np.array(low, np.int64)

rng = np.random.default_rng(316); random.seed(316)
chains = {a: chain(a) for a in (4, 5, 6, 7, 8)}
print("chains built", flush=True)
out = []
for bits, reps in ((100000, 4), (300000, 1)):
    for rep in range(reps):
        t0 = time.time(); y = (1 << bits) | random.getrandbits(bits - 1) | 1
        res, ks, low = orbit(y); n = len(ks); rec = dict(bits=bits, steps=n)
        # (1) total variation against pi_a, with i.i.d. calibration
        tv = {}
        for a in (5, 6, 7, 8):
            cls, idx, P, pi = chains[a]; r = res % 3 ** a
            emp = np.bincount(idx[r], minlength=len(cls)) / n
            d = 0.5 * np.abs(emp - pi).sum()
            sims = [0.5 * np.abs(np.bincount(rng.choice(len(cls), n, p=pi), minlength=len(cls)) / n - pi).sum() for _ in range(20)]
            tv[a] = dict(tv=d, iid_mean=float(np.mean(sims)), iid_sd=float(np.std(sims)))
        rec['tv'] = tv
        # (2) transitions mod 3^4 against the induced chain
        cls, idx, P, pi = chains[4]; r = idx[res % 81]
        C = np.zeros_like(P); np.add.at(C, (r[1:], r[:-1]), 1)
        def stat(C):
            E = P * C.sum(axis=0, keepdims=True); mask = E >= 5
            return float(((C - E) ** 2 / np.where(mask, E, 1))[mask].sum())
        chi2 = stat(C)
        # null distribution: the induced chain itself, simulated for the same length
        cum = np.cumsum(P, axis=0); sims = []
        for _ in range(40):
            s0 = rng.choice(P.shape[0], p=pi); u = rng.random(n); path = np.empty(n, np.int64); path[0] = s0
            for t in range(1, n): path[t] = np.searchsorted(cum[:, path[t - 1]], u[t] * cum[-1, path[t - 1]])
            Cs = np.zeros_like(P); np.add.at(Cs, (path[1:], path[:-1]), 1); sims.append(stat(Cs))
        sims = np.array(sims); dof = 0
        rec['transitions'] = dict(chi2=chi2, null_mean=float(sims.mean()), null_sd=float(sims.std()),
                                  p=float((np.sum(sims >= chi2) + 1) / (len(sims) + 1)))
        # (3) exponents
        freq = np.bincount(ks, minlength=12)[1:12] / n
        geo = 2.0 ** -np.arange(1, 12)
        kc = ks - ks.mean(); ac = [float((kc[:-L] * kc[L:]).mean() / kc.var()) for L in range(1, 6)]
        rec['k'] = dict(mean=float(ks.mean()), freq=freq.tolist(), geo=geo.tolist(), autocorr=ac, se=1 / math.sqrt(n),
                        z_mean=float((ks.mean() - 2) / math.sqrt(2 / n)))
        # (4) low binary digits
        cnt = np.bincount(low, minlength=4096)[1::2]; Eb = n / 2048
        chi2b = float(((cnt - Eb) ** 2 / Eb).sum())
        rec['low'] = dict(chi2=chi2b, dof=2047, p=float(stats.chi2.sf(chi2b, 2047)))
        out.append(rec)
        print(f"orbit from 2^{bits} ({n} steps, {time.time()-t0:.0f}s): TV/iid " +
              " ".join(f"a={a}: {tv[a]['tv']:.4f}/{tv[a]['iid_mean']:.4f}" for a in tv) +
              f"; transitions stat {chi2:.0f} vs chain null {rec['transitions']['null_mean']:.0f}+-{rec['transitions']['null_sd']:.0f}; mean k {ks.mean():.4f} (z={rec['k']['z_mean']:+.2f}); "
              f"autocorr {np.round(ac, 4)} (se {1/math.sqrt(n):.4f}); low digits chi2 {chi2b:.0f} on 2047 (p={rec['low']['p']:.3f})", flush=True)
json.dump(out, open(os.path.join(HERE, 'data', 'ch16_p3.json'), 'w'))
