import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
# Chapter 7, Section sec:ch7:experiments: the experiments with Switching MM, rerun for the book.
# Writes ch7_experiments.json (read by ch7_experiments_fig.py, ch7_noisefig.py and ch7_reconcile.py).
#
# Methods. Every method receives L and uses the step 1/L unless stated otherwise.
#   GD; SMM_det (Algorithm smmdet: Nesterov injection at odd t, plain gradient step at even t);
#   SMM_state (Algorithm smmstate with f_ref = f*, rho_thr = 0.356, eps = 1e-12);
#   fixed momentum heavy ball and Nesterov at beta in {0.5, log2(3) - 1, 2/3, 0.9};
#   tuned heavy ball (alpha = 4/(sqrt L + sqrt mu)^2, beta = ((sqrt L - sqrt mu)/(sqrt L + sqrt mu))^2);
#   tuned Nesterov (alpha = 1/L, beta = (sqrt L - sqrt mu)/(sqrt L + sqrt mu));
#   Adam (beta1 = 0.9, beta2 = 0.999, eps = 1e-8, step 1/L);
#   random injection (a Nesterov step with probability 1/2 at each t, a gradient step otherwise).
# Counts. The count of a run is the first t with f(x_t) - f* below the tolerance; a run that does not reach the
#   tolerance within the horizon T counts as T + 1. Means and standard errors are over R seeded replicates.
# Experiments.
#   1. Random quadratics f(x) = x^T A x / 2, d = 25, eigenvalues geometrically spaced in [1, kappa], random orthogonal
#      eigenbasis, x0 ~ N(0, I): iterations to f < 1e-4 for kappa in {2, 5, 10, 25, 50, 100}, T = 300, R = 20
#      (condition number sweep and frequency ablation). The runs are continued to T_FULL = 3000 steps, which
#      gives the counts without the cap at T + 1 and the number of runs that the cap censors.
#   2. Noise robustness: kappa = 25, gradients A x + sigma xi with xi ~ N(0, I) fresh at every call (the same draws
#      for every method within a replicate), sigma in {0, 0.05, 0.1, 0.2, 0.4}; final gap after T = 300 steps, R = 30.
#   3. l2 regularized logistic regression on the Breast Cancer Wisconsin data (standardized features, no intercept,
#      penalty 0.05), w0 ~ N(0, 0.1^2 I): iterations to a gap below 1e-4, T = 600, R = 20.
#   4. The same problem with minibatch gradients (batches drawn uniformly without replacement, the same batches for
#      every method within a replicate), |B| in {16, 64, 256}: iterations to a full gap below 1e-3, T = 600, R = 30.
import numpy as np, json, time
from scipy.special import expit
from sklearn.datasets import load_breast_cancer

B0 = np.log(1.5) / np.log(2)
BETAS = [0.5, B0, 2 / 3, 0.9]
BNAME = ['0.5', '0.585', '0.667', '0.9']
RHO_THR, EPS_STATE = 0.356, 1e-12
OUT = BOOK + '/scripts/ch7_experiments.json'
T0 = time.time()


def tuned(L, mu):
    sL, sm = np.sqrt(L), np.sqrt(mu)
    return {'hb': (4 / (sL + sm) ** 2, ((sL - sm) / (sL + sm)) ** 2), 'nag': (1 / L, (sL - sm) / (sL + sm))}


def col(v):
    v = np.asarray(v, float)
    return v.reshape(-1, 1) if v.ndim == 1 else v


def run(kind, beta, alpha, grad, fgap, x0, T, rng=None):
    """Run a method from the rows of x0 (one row per replicate). grad(x, t) is the (possibly noisy) gradient oracle,
    fgap(x) the objective gap of each row. Returns the gaps f(x_t) - f*, t = 0..T (shape (T+1, R)), and the number
    of injection steps of each run."""
    R = x0.shape[0]
    a, b = col(alpha), col(beta)
    x = x0.copy(); xp = x0.copy()
    gaps = np.empty((T + 1, R)); gaps[0] = fgap(x)
    ninj = np.zeros(R)
    if kind == 'adam':
        m = np.zeros_like(x); v = np.zeros_like(x)
    gprev = np.maximum(gaps[0], EPS_STATE)
    for t in range(T):
        if kind == 'hb':
            xn = x - a * grad(x, t) + b * (x - xp)
        elif kind == 'adam':
            g = grad(x, t); m = 0.9 * m + 0.1 * g; v = 0.999 * v + 0.001 * g * g
            xn = x - a * (m / (1 - 0.9 ** (t + 1))) / (np.sqrt(v / (1 - 0.999 ** (t + 1))) + 1e-8)
        else:
            if kind == 'gd':
                inj = np.zeros(R, bool)
            elif kind == 'nag':
                inj = np.ones(R, bool)
            elif kind == 'smm':
                inj = np.full(R, t % 2 == 1)
            elif kind == 'rand':
                inj = rng.random(R) < 0.5
            elif kind == 'state':
                gt = np.maximum(gaps[t], EPS_STATE)
                inj = np.full(R, t == 0) | (gt / gprev > RHO_THR)
                gprev = gt
            ninj += inj
            y = x + b * (x - xp) * inj[:, None]
            xn = y - a * grad(y, t)
        xp, x = x, xn
        gaps[t + 1] = fgap(x)
    return gaps, ninj


def count(gaps, tol):
    hit = gaps < tol
    return np.where(hit.any(0), hit.argmax(0), gaps.shape[0]).astype(float)


def ms(c):
    c = np.asarray(c, float)
    return [float(c.mean()), float(c.std(ddof=1) / np.sqrt(len(c)))]


def methods_fixed(L):
    out = {'GD': ('gd', 0.0, 1 / L), 'SMM': ('smm', B0, 1 / L)}
    for bb, nm in zip(BETAS, BNAME):
        out['HB ' + nm] = ('hb', bb, 1 / L)
        out['NAG ' + nm] = ('nag', bb, 1 / L)
    return out


res = {'setup': {'beta0': B0, 'rho_thr': RHO_THR}}

# ------------------------------------------------------------------ 1. random quadratics
def quad_batch(seed, R, d, kappa):
    rng = np.random.default_rng(seed)
    lam = np.geomspace(1.0, kappa, d); A = np.empty((R, d, d))
    for r in range(R):
        Q, Rm = np.linalg.qr(rng.standard_normal((d, d))); Q = Q * np.sign(np.diag(Rm))
        A[r] = (Q * lam) @ Q.T
    x0 = rng.standard_normal((R, d))
    return A, x0


def qgap(A):
    return lambda x: 0.5 * np.einsum('ri,rij,rj->r', x, A, x)


d, R, T, TOL, T_FULL = 25, 20, 300, 1e-4, 3000
res['ksweep'] = {'kappa': [2, 5, 10, 25, 50, 100], 'T': T, 'R': R, 'tol': TOL, 'T_full': T_FULL, 'counts': {},
                 'counts_full': {}, 'capped': {}}
for kappa in res['ksweep']['kappa']:
    A, x0 = quad_batch(100 + kappa, R, d, kappa)
    grad = lambda x, t, A=A: np.einsum('rij,rj->ri', A, x)
    meth = methods_fixed(float(kappa)); meth['random'] = ('rand', B0, 1 / kappa)
    for name, (kind, bb, al) in meth.items():
        g, _ = run(kind, bb, al, grad, qgap(A), x0, T_FULL, rng=np.random.default_rng(7000 + kappa))
        cf = count(g, TOL)
        res['ksweep']['counts'].setdefault(name, []).append(ms(np.minimum(cf, T + 1)))
        res['ksweep']['counts_full'].setdefault(name, []).append(ms(cf))
        res['ksweep']['capped'].setdefault(name, []).append(int((cf > T).sum()))
    print(f"kappa={kappa}: " + ", ".join(f"{k} {v[-1][0]:.1f}" for k, v in res['ksweep']['counts'].items()), flush=True)
    print(f"kappa={kappa} without the cap: " + ", ".join(
        f"{k} {v[-1][0]:.1f} ({res['ksweep']['capped'][k][-1]} of {R} capped)"
        for k, v in res['ksweep']['counts_full'].items()), flush=True)

# ------------------------------------------------------------------ 2. noise robustness, kappa = 25
R, T, kappa = 30, 300, 25.0
SIG = [0.0, 0.05, 0.10, 0.20, 0.40]
A, x0 = quad_batch(2025, R, d, kappa)
noise_meth = {'GD': ('gd', 0.0), 'SMM': ('smm', B0)}
for bb, nm in zip(BETAS, BNAME):
    noise_meth['NAG ' + nm] = ('nag', bb)
noise_meth['HB 0.585'] = ('hb', B0); noise_meth['HB 0.9'] = ('hb', 0.9)
res['noise'] = {'sigma': SIG, 'T': T, 'R': R, 'kappa': kappa, 'final': {}}
for j, s in enumerate(SIG):
    Xi = np.random.default_rng(3000 + j).standard_normal((T, R, d))
    grad = lambda x, t, s=s, Xi=Xi: np.einsum('rij,rj->ri', A, x) + s * Xi[t]
    for name, (kind, bb) in noise_meth.items():
        g, _ = run(kind, bb, 1 / kappa, grad, qgap(A), x0, T)
        res['noise']['final'].setdefault(name, []).append(ms(g[-1]))
    print(f"sigma={s}: " + ", ".join(f"{k} {v[-1][0]:.3e}" for k, v in res['noise']['final'].items()), flush=True)

# ------------------------------------------------------------------ 3. logistic regression
data = load_breast_cancer(); X = data.data; y = np.where(data.target == 1, 1.0, -1.0)
X = (X - X.mean(0)) / X.std(0); n, p = X.shape; LAM = 0.05


def lf(W):
    Z = (W @ X.T) * y
    return np.logaddexp(0.0, -Z).mean(1) + 0.5 * LAM * (W * W).sum(1)


def lgrad(W, t=None):
    Z = (W @ X.T) * y
    return -(expit(-Z) * y) @ X / n + LAM * W


w = np.zeros(p)
for _ in range(50):
    z = (X @ w) * y; s_ = expit(z)
    g = -(X.T @ (y * expit(-z))) / n + LAM * w
    H = (X.T * (s_ * (1 - s_))) @ X / n + LAM * np.eye(p)
    w = w - np.linalg.solve(H, g)
    if np.linalg.norm(g) < 1e-15: break
wstar = w; fstar = float(lf(wstar[None])[0])
z = (X @ wstar) * y; s_ = expit(z); Hs = (X.T * (s_ * (1 - s_))) @ X / n + LAM * np.eye(p)
L = float(np.linalg.eigvalsh(X.T @ X / n).max() / 4 + LAM); mu = LAM; kap = L / mu
tun = tuned(L, mu)
res['logistic'] = {'n': n, 'p': p, 'lam': LAM, 'L': L, 'mu': mu, 'kappa': kap, 'fstar': fstar,
                   'hess_min_over_L': float(np.linalg.eigvalsh(Hs).min() / L), 'T': 600, 'R': 20, 'tol': 1e-4,
                   'counts': {}, 'curves': {}, 'injections': {}}
R, T = 20, 600
W0 = 0.1 * np.random.default_rng(4000).standard_normal((R, p))
lgap = lambda W: lf(W) - fstar
meth = methods_fixed(L)
meth['HB tuned'] = ('hb',) + tun['hb'][::-1]
meth['NAG tuned'] = ('nag',) + tun['nag'][::-1]
meth['SMM state'] = ('state', B0, 1 / L)
meth['Adam'] = ('adam', 0.0, 1 / L)
for name, (kind, bb, al) in meth.items():
    g, ni = run(kind, bb, al, lgrad, lgap, W0, T)
    res['logistic']['counts'][name] = ms(count(g, 1e-4))
    res['logistic']['curves'][name] = np.maximum(g.mean(1), 1e-17).tolist()
    res['logistic']['injections'][name] = ms(ni)
    print(f"logistic {name}: {res['logistic']['counts'][name][0]:.1f} +- {res['logistic']['counts'][name][1]:.1f}"
          f" (injections {ni.mean():.1f} of {T})", flush=True)
print(f"logistic: L = {L:.4f}, kappa = {kap:.3f}, tuned NAG beta = {tun['nag'][1]:.4f}, f* = {fstar:.6f},"
      f" lambda_min(H*)/L = {res['logistic']['hess_min_over_L']:.5f}", flush=True)

# ------------------------------------------------------------------ 4. minibatch gradients
R, T = 30, 600
res['sgd'] = {'batch': [16, 64, 256], 'T': T, 'R': R, 'tol': 1e-3, 'counts': {}, 'curves': {}}
sgd_meth = {'GD': ('gd', 0.0, 1 / L), 'NAG tuned': ('nag', tun['nag'][1], 1 / L),
            'NAG 0.585': ('nag', B0, 1 / L), 'NAG 0.9': ('nag', 0.9, 1 / L), 'SMM': ('smm', B0, 1 / L)}
W0 = 0.1 * np.random.default_rng(5000).standard_normal((R, p))
for Bsz in res['sgd']['batch']:
    rng = np.random.default_rng(6000 + Bsz)
    IDX = np.stack([rng.random((R, n)).argsort(1)[:, :Bsz] for _ in range(T)])

    def sgrad(W, t, IDX=IDX, Bsz=Bsz):
        Xb = X[IDX[t]]; yb = y[IDX[t]]
        Z = np.einsum('rbp,rp->rb', Xb, W) * yb
        return -np.einsum('rbp,rb->rp', Xb, expit(-Z) * yb) / Bsz + LAM * W
    for name, (kind, bb, al) in sgd_meth.items():
        g, _ = run(kind, bb, al, sgrad, lgap, W0, T)
        c = count(g, 1e-3)
        res['sgd']['counts'].setdefault(name, []).append(ms(c) + [int((c > T).sum())])
        res['sgd']['curves'].setdefault(name, []).append(g.mean(1).tolist())
        res['sgd'].setdefault('plateau', {}).setdefault(name, []).append(float(g[-200:].mean()))
    print(f"|B|={Bsz}: " + ", ".join(f"{k} {v[-1][0]:.1f}+-{v[-1][1]:.1f} ({v[-1][2]} capped)"
                                      for k, v in res['sgd']['counts'].items()), flush=True)
    print(f"|B|={Bsz} plateaus (mean gap over the last 200 iterations): " +
          ", ".join(f"{k} {v[-1]:.2e}" for k, v in res['sgd']['plateau'].items()), flush=True)

json.dump(res, open(OUT, 'w'))
print(f"done in {time.time() - T0:.1f} s", flush=True)
