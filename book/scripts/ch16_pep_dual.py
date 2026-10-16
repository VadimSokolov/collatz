"""P11: dual multipliers of the convex performance estimation problem for Switching MM (cyclic GD with steps
h = 1 + beta, 1).  Points: star (index -1 shown as '*'), p_0 .. p_N.  Inequalities
  f_i >= f_j + <g_j, x_i - x_j> + |g_i - g_j|^2 / 2      (L = 1)
with multipliers lambda_{ij}.  Prints the multipliers above a threshold, to guess the closed form."""
import numpy as np, cvxpy as cp, math, sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ch5d_pep import smm, cyclic_gd

def pep_dual(H):
    N = H.shape[0] - 1; dim = 1 + (N + 1); E = np.eye(dim)
    P = [E[0] - sum(H[k, j] * E[1 + j] for j in range(k)) for k in range(N + 1)]
    Gr = [E[1 + k] for k in range(N + 1)]
    G = cp.Variable((dim, dim), PSD=True); F = cp.Variable(N + 1)
    pts = [(np.zeros(dim), np.zeros(dim), 0.0)] + [(P[k], Gr[k], F[k]) for k in range(N + 1)]
    names = ['*'] + [str(k) for k in range(N + 1)]
    ip = lambda a, b: a @ G @ b
    cons = [E[0] @ G @ E[0] <= 1]; idx = []
    for i in range(len(pts)):
        for j in range(len(pts)):
            if i == j: continue
            xi, gi, fi = pts[i]; xj, gj, fj = pts[j]
            cons.append(fi >= fj + ip(gj, xi - xj) + 0.5 * ip(gi - gj, gi - gj)); idx.append((names[i], names[j]))
    prob = cp.Problem(cp.Maximize(F[N]), cons)
    prob.solve(solver='CLARABEL', tol_gap_abs=1e-10, tol_gap_rel=1e-10, tol_feas=1e-10)
    lam = {idx[t]: cons[1 + t].dual_value for t in range(len(idx))}
    return prob.value, cons[0].dual_value, lam

if __name__ == '__main__':
    out = {}
    for beta in (0.0, 0.5, math.log(1.5) / math.log(2)):
        for N in (2, 4, 6):
            val, tau, lam = pep_dual(smm(N, beta))
            big = {k: float(v) for k, v in lam.items() if abs(v) > 1e-6}
            print(f"beta={beta:.4f} N={N}: value {val:.8f} (formula {1/((4+2*beta)*N+2):.8f}), dual of radius {float(tau):.8f}")
            for k in sorted(big, key=lambda t: (t[0] != '*', t)):
                print(f"   lambda[{k[0]},{k[1]}] = {big[k]:.8f}")
            out[f"{beta:.4f}_{N}"] = dict(value=val, tau=float(tau), lam={f"{a},{b}": v for (a, b), v in big.items()})
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'ch16_pep_dual.json'), 'w'))
