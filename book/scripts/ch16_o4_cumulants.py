"""Open problem O4: the third cumulant of total stopping times.

(1) Census (exact, n < 2^30, from data/ch3_decomp.npz): second and third cumulants of tau over each
    dyadic block, against the Bernoulli model.  Least squares fit kappa_j(k) = A_j E[log n] - D_j over
    k >= 20, and comparison of the slopes A_j with the model values sigma^2/|mu|^3 and 3 sigma^4/|mu|^5.
(2) Split at Y = 2^14.  For n uniform in a block, real orbits are split at the first time below Y:
    tau = tau_top + tau_bot, tau_bot = tau(y) read from the census table.  Model walks (steps
    log(3/2), -log 2 with probability 1/2, started at log n, stopped below 0) are split in the same way.
    Cumulants of each part and of the cross terms, for k = 20, 24, 28.
Output data/ch16_o4.json."""
import numpy as np, numba as nb, json, math, os
HERE = os.path.dirname(os.path.abspath(__file__))
L2, L3, L32 = math.log(2), math.log(3), math.log(1.5)
MU = abs(math.log(math.sqrt(3) / 2)); S2 = (L3 / 2) ** 2
A2, A3 = S2 / MU ** 3, 3 * S2 ** 2 / MU ** 5

rows = np.load(os.path.join(HERE, 'data', 'ch3_decomp.npz'))['rows']
ks = rows[:, 0].astype(int); eu = rows[:, 2]; k3 = rows[:, 7]; k2 = rows[:, 9]
sel = ks >= 20
fit = {}
for name, y, A in (('k2', k2, A2), ('k3', k3, A3)):
    X = np.vstack([eu[sel], -np.ones(sel.sum())]).T
    (a, d), *_ = np.linalg.lstsq(X, y[sel], rcond=None)
    resid = y[sel] - X @ np.array([a, d])
    fit[name] = dict(slope=a, model_slope=A, ratio=a / A, D=d, max_resid=float(np.abs(resid).max()))
    print(f"{name}: fitted slope {a:.2f} vs model {A:.2f} (ratio {a/A:.4f}); deficit D = {d:.1f}; max residual {np.abs(resid).max():.2f}")
# predicted skewness from the fit, and observed
skew_obs = k3 / k2 ** 1.5
pred = (fit['k3']['slope'] * eu - fit['k3']['D']) / (fit['k2']['slope'] * eu - fit['k2']['D']) ** 1.5
model_skew = (A3 * eu) / (A2 * eu) ** 1.5
for k, so, sp, sm in zip(ks, skew_obs, pred, model_skew):
    if k >= 12: print(f"  k={k}: skew observed {so:.4f}, fit {sp:.4f}, model {sm:.4f}")
future = {}
for kk in (40, 50, 60, 100):
    u = (kk + 0.5) * L2
    future[kk] = dict(fit=(fit['k3']['slope'] * u - fit['k3']['D']) / (fit['k2']['slope'] * u - fit['k2']['D']) ** 1.5,
                      model=A3 * u / (A2 * u) ** 1.5)
    print(f"  extrapolated k={kk}: skew {future[kk]['fit']:.4f} vs model {future[kk]['model']:.4f}")

tau = np.load(os.path.join(HERE, 'data', 'ch3_tau_small.npy')).astype(np.int64)
@nb.njit(cache=True)
def real_split(ns, Y, tau):
    m = ns.shape[0]; top = np.zeros(m, np.int64); bot = np.zeros(m, np.int64)
    for i in range(m):
        x = ns[i]; s = 0
        while x >= Y:
            if x & 1: x = (3 * x + 1) >> 1
            else: x >>= 1
            s += 1
        top[i] = s; bot[i] = tau[x]
    return top, bot
@nb.njit(cache=True)
def model_split(us, logY, seed):
    np.random.seed(seed)
    m = us.shape[0]; top = np.zeros(m, np.int64); bot = np.zeros(m, np.int64)
    a, b = math.log(1.5), math.log(2.0)
    for i in range(m):
        u = us[i]; s = 0
        while u >= logY:
            u = u + a if np.random.random() < 0.5 else u - b
            s += 1
        t = 0
        while u >= 0.0:
            u = u + a if np.random.random() < 0.5 else u - b
            t += 1
        top[i] = s; bot[i] = t
    return top, bot
def cum(x):
    x = x.astype(float); c = x - x.mean(); return c.var(), (c ** 3).mean()
def cross3(a, b):
    a = a - a.mean(); b = b - b.mean(); return float((a * a * b).mean()), float((a * b * b).mean())
rng = np.random.default_rng(16); Y = 1 << 14; split = {}
for k in (20, 24, 28):
    ns = rng.integers(1 << k, 1 << (k + 1), 1_000_000).astype(np.int64)
    rt, rb = real_split(ns, Y, tau)
    mt, mb = model_split(np.log(ns.astype(float)), math.log(Y), 100 + k)
    out = {}
    for nm, t, b in (('real', rt, rb), ('model', mt, mb)):
        v_t, c_t = cum(t); v_b, c_b = cum(b); v, c = cum(t + b)
        x1, x2 = cross3(t, b)
        out[nm] = dict(var_top=v_t, k3_top=c_t, var_bot=v_b, k3_bot=c_b, var=v, k3=c, cov=float(np.cov(t, b)[0, 1]),
                       cross=3 * (x1 + x2), mean_bot=float(b.mean()))
    split[k] = out
    R, M = out['real'], out['model']
    print(f"k={k}: deficit var total {M['var']-R['var']:.1f} = top {M['var_top']-R['var_top']:.1f} + bottom {M['var_bot']-R['var_bot']:.1f} + 2cov {2*(M['cov']-R['cov']):.1f}")
    print(f"      deficit k3  total {M['k3']-R['k3']:.0f} = top {M['k3_top']-R['k3_top']:.0f} + bottom {M['k3_bot']-R['k3_bot']:.0f} + cross {M['cross']-R['cross']:.0f}")
    print(f"      real bottom: mean {R['mean_bot']:.2f} var {R['var_bot']:.1f} k3 {R['k3_bot']:.0f}; model bottom: mean {M['mean_bot']:.2f} var {M['var_bot']:.1f} k3 {M['k3_bot']:.0f}")
json.dump(dict(fit=fit, skew=dict(k=ks.tolist(), obs=skew_obs.tolist(), fit=pred.tolist(), model=model_skew.tolist()),
               future={str(k): v for k, v in future.items()}, split={str(k): v for k, v in split.items()}),
          open(os.path.join(HERE, 'data', 'ch16_o4.json'), 'w'), default=float)
