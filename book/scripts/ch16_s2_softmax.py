"""S2 and P5: noise floors of Switching MM on a real stochastic problem.
L2 regularized multinomial logistic regression on the scikit learn digits data (1797 x 64, ten classes,
pixels scaled to [0,1], a bias feature), lambda = 0.01.  The minimizer w_* is computed by L-BFGS; H is the
Hessian at w_*, L its largest eigenvalue, S the covariance of per sample gradients at w_*.
Prediction: the linearized dynamics around w_* with minibatch noise of covariance S/b give stationary
covariances by discrete Lyapunov equations, and the excess loss tr(H C)/2, for gradient descent, Switching MM
(beta = log2 3 - 1), every step Nesterov at beta = 0.585 and 0.9, and heavy ball at beta = 0.9, all with
step 1/L.  Simulation: each method run on the true loss from w_* with minibatches of size b (with
replacement), excess loss averaged over a long horizon and several seeds.  Output data/ch16_s2.json."""
import numpy as np, json, math, os, time
from scipy.optimize import minimize
from scipy.linalg import solve_discrete_lyapunov
from sklearn.datasets import load_digits
HERE = os.path.dirname(os.path.abspath(__file__))
X, y = load_digits(return_X_y=True); X = np.hstack([X / 16.0, np.ones((X.shape[0], 1))]); n, p = X.shape; C = 10
Y = np.eye(C)[y]; LAM = 0.01; BETA = math.log(1.5) / math.log(2)
def softmax(Z):
    Z = Z - Z.max(axis=1, keepdims=True); E = np.exp(Z); return E / E.sum(axis=1, keepdims=True)
def loss(w, idx=None):
    W = w.reshape(p, C); Xb, Yb = (X, Y) if idx is None else (X[idx], Y[idx])
    Z = Xb @ W; Z = Z - Z.max(axis=1, keepdims=True); lse = np.log(np.exp(Z).sum(axis=1))
    return float((lse - (Z * Yb).sum(axis=1)).mean() + 0.5 * LAM * (w @ w))
def grad(w, idx=None):
    W = w.reshape(p, C); Xb, Yb = (X, Y) if idx is None else (X[idx], Y[idx])
    P = softmax(Xb @ W); return ((Xb.T @ (P - Yb)) / Xb.shape[0]).ravel() + LAM * w
r = minimize(lambda w: (loss(w), grad(w)), np.zeros(p * C), jac=True, method='L-BFGS-B', options=dict(maxiter=5000, gtol=1e-12))
ws = r.x; fs = loss(ws); print(f"f* = {fs:.8f}, |grad| = {np.linalg.norm(grad(ws)):.2e}", flush=True)
# Hessian at w_*
P = softmax(X @ ws.reshape(p, C)); d = p * C; H = np.zeros((d, d))
for i in range(n):
    Pi = np.diag(P[i]) - np.outer(P[i], P[i]); H += np.kron(np.outer(X[i], X[i]), Pi)
H = H / n + LAM * np.eye(d); H = (H + H.T) / 2
ev = np.linalg.eigvalsh(H); L = ev.max(); mu = ev.min(); print(f"L = {L:.4f}, mu = {mu:.4f}, kappa = {L/mu:.1f}", flush=True)
# per sample gradient covariance
G = np.einsum('ij,ik->ijk', X, P - Y).reshape(n, d) + LAM * ws; Gc = G - G.mean(axis=0); S = Gc.T @ Gc / n
alpha = 1.0 / L; I = np.eye(d); Z0 = np.zeros((d, d)); M = I - alpha * H
def blocks(beta):
    A = np.block([[(1 + beta) * M, -beta * M], [I, Z0]]); return A
def hb(beta):
    return np.block([[M + beta * I, -beta * I], [I, Z0]])
Bn = np.vstack([-alpha * I, Z0])
def predict(b):
    Q = S / b; out = {}
    Cg = solve_discrete_lyapunov(M, alpha ** 2 * Q); out['gd'] = 0.5 * np.trace(H @ Cg)
    for nm, A in (('nag585', blocks(BETA)), ('nag9', blocks(0.9)), ('hb9', hb(0.9))):
        Cz = solve_discrete_lyapunov(A, Bn @ Q @ Bn.T); out[nm] = 0.5 * np.trace(H @ Cz[:d, :d])
    Ae = np.block([[M, Z0], [I, Z0]]); Ao = blocks(BETA)
    A2 = Ao @ Ae; Q2 = Ao @ Bn @ Q @ Bn.T @ Ao.T + Bn @ Q @ Bn.T
    Ce = solve_discrete_lyapunov(A2, Q2); Co = Ae @ Ce @ Ae.T + Bn @ Q @ Bn.T
    out['smm'] = 0.25 * (np.trace(H @ Ce[:d, :d]) + np.trace(H @ Co[:d, :d]))
    return out
def simulate(method, b, steps, seed, burn=3000, every=10):
    rng = np.random.default_rng(seed); w = ws.copy(); wp = ws.copy(); acc = []; t0 = 0
    for t in range(burn + steps):
        idx = rng.integers(0, n, b)
        if method == 'gd' or (method == 'smm' and t % 2 == 0):
            wn = w - alpha * grad(w, idx)
        elif method in ('nag585', 'nag9') or method == 'smm':
            beta = 0.9 if method == 'nag9' else BETA
            yv = w + beta * (w - wp); wn = yv - alpha * grad(yv, idx)
        elif method == 'hb9':
            wn = w - alpha * grad(w, idx) + 0.9 * (w - wp)
        wp, w = w, wn
        if t >= burn and t % every == 0: acc.append(loss(w) - fs)
    return float(np.mean(acc))
res = {}
for b in (16, 64):
    pr = predict(b); sims = {}
    for m in ('gd', 'smm', 'nag585', 'nag9', 'hb9'):
        vals = [simulate(m, b, 20000, 1000 * b + s) for s in range(6)]
        sims[m] = (float(np.mean(vals)), float(np.std(vals) / math.sqrt(len(vals))))
    res[b] = dict(pred=pr, sim=sims)
    print(f"b={b}:", flush=True)
    for m in pr:
        print(f"   {m:7s}: predicted {pr[m]:.3e}  simulated {sims[m][0]:.3e} +- {sims[m][1]:.1e}  ratio to GD: pred {pr[m]/pr['gd']:.3f} sim {sims[m][0]/sims['gd'][0]:.3f}", flush=True)
json.dump(dict(L=L, mu=mu, kappa=L / mu, fstar=fs, res={str(k): v for k, v in res.items()}),
          open(os.path.join(HERE, 'data', 'ch16_s2.json'), 'w'), default=float)
