"""P11: sparse certificate.  Keep only the inequalities (i, i+1), (i+1, i) between consecutive points and
(*, i) between the minimizer and each point.  If the relaxed program has the same value as the full one,
a certificate with these multipliers exists; print them."""
import numpy as np, cvxpy as cp, math, sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ch5d_pep import smm, pep
def pep_sparse(H, pairs):
    N = H.shape[0] - 1; dim = 1 + (N + 1); E = np.eye(dim)
    P = [E[0] - sum(H[k, j] * E[1 + j] for j in range(k)) for k in range(N + 1)]
    Gr = [E[1 + k] for k in range(N + 1)]
    G = cp.Variable((dim, dim), PSD=True); F = cp.Variable(N + 1)
    pt = {'*': (np.zeros(dim), np.zeros(dim), 0.0)}
    for k in range(N + 1): pt[k] = (P[k], Gr[k], F[k])
    ip = lambda a, b: a @ G @ b
    cons = [E[0] @ G @ E[0] <= 1]
    for (i, j) in pairs:
        xi, gi, fi = pt[i]; xj, gj, fj = pt[j]
        cons.append(fi >= fj + ip(gj, xi - xj) + 0.5 * ip(gi - gj, gi - gj))
    prob = cp.Problem(cp.Maximize(F[N]), cons); prob.solve(solver='CLARABEL')
    return prob.value, float(cons[0].dual_value), {pairs[t]: float(cons[1 + t].dual_value) for t in range(len(pairs))}
B = math.log(1.5) / math.log(2)
for beta in (0.0, 0.5, B):
    for N in (2, 4, 6, 8):
        H = smm(N, beta); full = 1 / ((4 + 2 * beta) * N + 2)
        sets = {'consecutive + star': [('*', k) for k in range(N + 1)] + [(k, k + 1) for k in range(N)] + [(k + 1, k) for k in range(N)],
                'forward + star': [('*', k) for k in range(N + 1)] + [(k, k + 1) for k in range(N)],
                'backward + star': [('*', k) for k in range(N + 1)] + [(k + 1, k) for k in range(N)]}
        for nm, pairs in sets.items():
            v, tau, lam = pep_sparse(H, pairs)
            ok = abs(v - full) / full < 1e-6
            print(f"beta={beta:.3f} N={N} {nm:20s}: value {v:.8f} vs formula {full:.8f} {'EXACT' if ok else ''}")
            if ok and nm != 'consecutive + star':
                print('    ', ' '.join(f"{a},{b}:{l:.6f}" for (a, b), l in lam.items() if abs(l) > 1e-7))
