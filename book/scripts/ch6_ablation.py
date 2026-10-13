# Chapter 6: the frequency ablation of Chapter 7 at beta = log2(3) - 1, run to convergence (censoring check).
# Random quadratics as in Chapter 7: d = 25, eigenvalues geometrically spaced in [1, kappa], step alpha = 1/L,
# x0 drawn from N(0, I) (isotropic, so the random orthogonal eigenbasis is immaterial), stop at an objective gap
# below 1e-4. Prints the mean iteration counts run to convergence and the same means censored at the horizon
# T = 300 of Chapter 7, together with the share of runs that had not converged by T.
import numpy as np
b = np.log(1.5) / np.log(2); d = 25; R = 200; T = 300; tol = 1e-4
rng = np.random.default_rng(2026)

def iterations(lam, e0, kind):
    e = e0.copy(); ep = e0.copy(); al = 1 / lam.max()
    for t in range(100000):
        if 0.5 * np.sum(lam * e ** 2) < tol:
            return t
        if kind == 'every': inject = True
        elif kind == 'parity': inject = (t % 2 == 1)
        elif kind == 'random': inject = rng.random() < 0.5
        else: inject = False
        y = e + b * (e - ep) if inject else e
        ep, e = e, y - al * lam * y
    return 100000

for kappa in [2, 10, 100]:
    lam = np.geomspace(1, kappa, d)
    counts = {k: [] for k in ['every', 'parity', 'random', 'never']}
    for r in range(R):
        e0 = rng.standard_normal(d)
        for k in counts:
            counts[k].append(iterations(lam, e0, k))
    print('kappa = %d' % kappa)
    for k, c in counts.items():
        c = np.array(c)
        print('  %-6s mean %6.1f   censored at T: %6.1f   share not converged by T: %.2f'
              % (k, c.mean(), np.minimum(c, T).mean(), np.mean(c > T)))
    print('  ratio parity/every: %.3f run to convergence, %.3f censored'
          % (np.mean(counts['parity']) / np.mean(counts['every']),
             np.minimum(counts['parity'], T).mean() / np.minimum(counts['every'], T).mean()))
