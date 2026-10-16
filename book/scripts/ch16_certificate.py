"""P11: closed form certificate for the convex worst case of Switching MM.
Cyclic gradient method with steps s_0 = h, s_1 = 1, s_2 = h, ... (h = 1 + beta), N even, L = 1, x_* = 0.
Multipliers: W_k = 2 (s_0 + ... + s_k), D = (4 + 2 beta) N + 2,
   lam_k = lam_{k,k+1} = W_k / (D - W_k + 2 beta [k even]),
   mu_0 = lam_0, mu_k = lam_k - lam_{k-1} (0 < k < N), mu_N = 1 - lam_{N-1}   (minimizer multipliers).
The weighted sum of the inequalities
   (k,k+1):  f_{k+1} - f_k + <g_{k+1}, x_k - x_{k+1}> + |g_k - g_{k+1}|^2/2 <= 0
   (*,k)  :  f_k - f_*   + <g_k, x_* - x_k>        + |g_k|^2/2           <= 0
equals f_N - f_* + Q(x_0, g).  The bound f_N - f_* <= |x_0|^2 / D holds if S = |x_0|^2/D + Q is PSD.
We compute S symbolically in beta for N = 2, 4, 6 and numerically for larger N, and locate where it
stops being PSD."""
import sympy as sp, numpy as np, math, json, os

def certificate_matrix(N, beta, symbolic=True):
    num = (lambda v: v) if symbolic else float
    h = 1 + beta
    s = [h if k % 2 == 0 else sp.Integer(1) if symbolic else 1.0 for k in range(N)]
    W = []; acc = 0
    for k in range(N): acc = acc + 2 * s[k]; W.append(acc)
    D = (4 + 2 * beta) * N + 2
    lam = [W[k] / (D - W[k] + (2 * beta if k % 2 == 0 else 0)) for k in range(N)]
    mu = [lam[0]] + [lam[k] - lam[k - 1] for k in range(1, N)] + [1 - lam[N - 1]]
    dim = N + 2                                  # basis x0, g_0..g_N
    Z = sp.zeros(dim, dim) if symbolic else np.zeros((dim, dim))
    def vec(d):
        v = sp.zeros(dim, 1) if symbolic else np.zeros(dim)
        for j, c in d.items(): v[j] = v[j] + c
        return v
    X = [vec({0: 1})]
    for k in range(N): X.append(X[-1] - s[k] * vec({1 + k: 1}))
    Gv = [vec({1 + k: 1}) for k in range(N + 1)]
    outer = (lambda a, b: a * b.T) if symbolic else (lambda a, b: np.outer(a, b))
    sym = lambda M: (M + M.T) / 2
    Q = Z.copy() if symbolic else Z.copy()
    for k in range(N):
        Q = Q + lam[k] * (sym(outer(Gv[k + 1], X[k] - X[k + 1])) + outer(Gv[k] - Gv[k + 1], Gv[k] - Gv[k + 1]) / 2)
    for k in range(N + 1):
        Q = Q + mu[k] * (sym(outer(Gv[k], -X[k])) + outer(Gv[k], Gv[k]) / 2)
    S = Q + outer(X[0], X[0]) / D
    return S, lam, mu, D

if __name__ == '__main__':
    b = sp.symbols('beta', positive=True)
    res = {}
    for N in (2, 4):
        S, lam, mu, D = certificate_matrix(N, b)
        S = sp.simplify(S)
        # S should be singular (the bound is attained); factor the leading principal minors
        minors = [sp.factor(sp.simplify(S[:j, :j].det())) for j in range(1, S.shape[0] + 1)]
        print(f"N={N}: rank of S at beta=0.585: {np.linalg.matrix_rank(np.array(S.subs(b, 0.585), dtype=float), tol=1e-9)}")
        for j, m in enumerate(minors, 1):
            print(f"   leading minor {j}: {m}")
        ev = np.linalg.eigvalsh(np.array(S.subs(b, sp.Rational(585, 1000)), dtype=float))
        print("   eigenvalues at beta=0.585:", np.round(ev, 10))
        res[N] = [str(m) for m in minors]
    # numerical: smallest eigenvalue of S as a function of beta, for N up to 20
    num = {}
    for N in (2, 4, 6, 8, 10, 14, 20):
        lo, hi = 0.5, 0.99
        def mineig(beta):
            S, *_ = certificate_matrix(N, beta, symbolic=False); return np.linalg.eigvalsh(S).min()
        betas = np.linspace(0.0, 0.98, 99); m = np.array([mineig(x) for x in betas])
        ok = betas[m > -1e-10]
        # threshold: largest beta in [0, 0.98] up to which S is PSD
        thr = None
        for x, v in zip(betas, m):
            if v < -1e-10: thr = x; break
        if thr is not None:
            a, c = thr - 0.01, thr
            for _ in range(50):
                mid = (a + c) / 2
                if mineig(mid) < -1e-12: c = mid
                else: a = mid
            thr = (a + c) / 2
        num[N] = thr; print(f"N={N}: S is PSD for beta up to {thr}")
    json.dump(dict(minors=res, thresholds={str(k): v for k, v in num.items()}),
              open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'ch16_certificate.json'), 'w'))
