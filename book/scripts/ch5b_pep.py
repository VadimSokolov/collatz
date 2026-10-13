"""Chapter 5b: performance estimation (PEP) for two steps of SMM_det on F_{mu,L}.
Even subsequence: y = x0 - (1+beta)/L g0, x2 = y - (1/L) g_y  (Lemma ch6:ystep).
Interpolation conditions of Taylor, Hendrickx and Glineur for smooth strongly convex functions.
Normalization L = 1, mu = 1/kappa."""
import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
import numpy as np, cvxpy as cp

def interp(fi, fj, xi, xj, gi, gj, G, mu, L):
    # f_i >= f_j + <g_j, x_i - x_j> + 1/(2(1-mu/L)) (|g_i-g_j|^2/L + mu|x_i-x_j|^2 - 2 mu/L <g_j-g_i, x_j-x_i>)
    ip = lambda a, b: a @ G @ b
    dx = xi - xj; dg = gi - gj
    return fi >= fj + ip(gj, dx) + 1/(2*(1-mu/L))*(ip(dg, dg)/L + mu*ip(dx, dx) - 2*mu/L*ip(-dg, -dx))

def pep(kappa, steps, objective='dist'):
    """steps: list of step lengths h (times 1/L); returns worst case ratio"""
    L, mu = 1.0, 1.0/kappa
    n = len(steps) + 1            # gradients g at x0, x1, ..., x_n-1 (and final point if fval)
    dim = 1 + n                   # basis: e0, g_0 .. g_{n-1}  (g at final point is basis vector n)
    G = cp.Variable((dim, dim), PSD=True)
    F = cp.Variable(n)            # f values at x0..x_{n-1}... plus final
    E = np.eye(dim)
    xs = [E[0]]; gs = [E[1 + i] for i in range(n)]
    for i, h in enumerate(steps):
        xs.append(xs[-1] - h*gs[i])
    # points: star (x=0,g=0,f=0) and x_0..x_K with gradients g_0..g_K (K = len(steps))
    pts = [(np.zeros(dim), np.zeros(dim), 0.0)] + [(xs[i], gs[i], F[i]) for i in range(n)]
    cons = []
    for i in range(len(pts)):
        for j in range(len(pts)):
            if i != j:
                cons.append(interp(pts[i][2], pts[j][2], pts[i][0], pts[j][0], pts[i][1], pts[j][1], G, mu, L))
    if objective == 'dist':
        cons.append(E[0] @ G @ E[0] <= 1)
        xK = xs[-1]; obj = xK @ G @ xK
    else:
        cons.append(F[0] <= 1)
        obj = F[n-1]
    prob = cp.Problem(cp.Maximize(obj), cons)
    prob.solve(solver='CLARABEL')
    return prob.value

if __name__ == '__main__':
    # sanity: one gradient step, distance: max(|1-h mu|,|1-h L|)^2
    for kappa, h in [(10, 1.0), (10, 1.5), (4, 1.8)]:
        print('GD check', kappa, h, pep(kappa, [h]), max(abs(1-h/kappa), abs(1-h))**2)
    print('GD one step fval', pep(10, [1.0], 'f'), (1-1/10)**2, ' h=1.585:', pep(10, [1.585], 'f'), max((1-1.585/10)**2, 0.585**2))
    beta = np.log(1.5)/np.log(2)
    rows = []
    for b in [0.3, beta, 0.9]:
        for kappa in [2, 3, 5, 10, 20, 50, 100]:
            q = (1-1/kappa)*max(b, 1-(1+b)/kappa)
            qf = (1-1/kappa)*(1-(1-b*b)/kappa)
            # two steps of SMM on even subsequence: long step then short step
            pd = pep(kappa, [1+b, 0.0, 1.0][:1] + [1.0])     # x0 -> y (h=1+b) -> x2 (h=1)
            pf = pep(kappa, [1+b, 1.0], 'f')
            pg = pep(kappa, [1.0, 1.0], 'f')
            rows.append((b, kappa, q*q, pd, qf, pf, (1-1/kappa)**2, pg))
            print(f"beta={b:.3f} kappa={kappa:4d}  q^2={q*q:.6f} PEPdist={pd:.6f} | qf={qf:.6f} PEPf={pf:.6f} q^2={q*q:.6f} | GD2 PL={(1-1/kappa)**2:.6f} PEP_GD2_f={pg:.6f}", flush=True)
    np.save(BOOK + '/scripts/ch5b_pep.npy', np.array(rows))
