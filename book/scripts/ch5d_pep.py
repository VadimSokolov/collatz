"""Performance estimation on smooth convex functions (mu = 0, L = 1) for first order methods
written in evaluation point form: p_0 = x_0, p_k = x_0 - sum_{j<k} H[k,j] g_j, gradients g_k at p_k,
output point p_N.  Worst case of f(p_N) - f* subject to |x_0 - x*| <= 1.
Interpolation conditions of Taylor, Hendrickx and Glineur for F_{0,L}:
  f_i >= f_j + <g_j, x_i - x_j> + |g_i - g_j|^2 / (2L).
Methods: GD, cyclic GD (Switching MM is cyclic GD with steps 1+beta, 1), constant momentum Nesterov,
Nesterov's method of similar triangles (FGM), and switched similar triangles (accelerated on odd
iterations, plain gradient steps on even ones)."""
import numpy as np, cvxpy as cp, json, math, sys, os

BETA = math.log(1.5) / math.log(2)

def pep(H):
    N = H.shape[0] - 1                  # points p_0..p_N, gradients g_0..g_N
    dim = 1 + (N + 1)                   # basis x0 (relative to x*=0), g_0..g_N
    E = np.eye(dim)
    P = [E[0] - sum(H[k, j] * E[1 + j] for j in range(k)) for k in range(N + 1)]
    Gr = [E[1 + k] for k in range(N + 1)]
    G = cp.Variable((dim, dim), PSD=True); F = cp.Variable(N + 1)
    pts = [(np.zeros(dim), np.zeros(dim), 0.0)] + [(P[k], Gr[k], F[k]) for k in range(N + 1)]
    ip = lambda a, b: a @ G @ b
    cons = [E[0] @ G @ E[0] <= 1]
    for i in range(len(pts)):
        for j in range(len(pts)):
            if i == j: continue
            xi, gi, fi = pts[i]; xj, gj, fj = pts[j]
            cons.append(fi >= fj + ip(gj, xi - xj) + 0.5 * ip(gi - gj, gi - gj))
    prob = cp.Problem(cp.Maximize(F[N]), cons)
    prob.solve(solver='CLARABEL')
    return prob.value

def H_from_points(points):
    """points: list of length N+1, each a dict {j: coeff} meaning p_k = x0 - sum coeff_j g_j."""
    N = len(points) - 1; H = np.zeros((N + 1, N + 1))
    for k, d in enumerate(points):
        for j, c in d.items(): H[k, j] = c
    return H

def cyclic_gd(steps):
    pts = [dict()]
    for h in steps:
        d = dict(pts[-1]); k = len(pts) - 1; d[k] = d.get(k, 0) + h; pts.append(d)
    return H_from_points(pts)

def smm(N, beta=BETA):
    return cyclic_gd([(1 + beta) if (t % 2 == 0) else 1.0 for t in range(N)])

def lin(*terms):
    out = {}
    for c, d in terms:
        for j, v in d.items(): out[j] = out.get(j, 0) + c * v
    return out

def nag_const(N, beta):
    # y_k = x_k + beta (x_k - x_{k-1}), x_{k+1} = y_k - g(y_k); evaluation points y_0..y_{N-1}, output x_N
    x_prev = dict(); x = dict(); pts = []
    for k in range(N):
        y = lin((1 + beta, x), (-beta, x_prev)); pts.append(y)
        xn = lin((1, y), (1, {k: 1.0})); x_prev, x = x, xn
    pts.append(x)
    return H_from_points(pts)

def similar_triangles(N, accel_mask):
    """accel_mask[t] True: similar triangles step (gradient at y); False: plain step at x."""
    A = 0.0; x = dict(); z = dict(); pts = []
    for t in range(N):
        if accel_mask[t]:
            a = (1 + math.sqrt(1 + 4 * A)) / 2; An = A + a
            y = lin((A / An, x), (a / An, z)); pts.append(y)
            z = lin((1, z), (a, {t: 1.0}))            # z+ = z - a g(y)
            x = lin((1, y), (1, {t: 1.0}))            # x+ = y - g(y)
            A = An
        else:
            pts.append(dict(x)); x = lin((1, x), (1, {t: 1.0}))
    pts.append(x)
    return H_from_points(pts), A

if __name__ == '__main__':
    Ns = list(range(1, 21))
    res = {k: [] for k in ('gd', 'smm', 'nag585', 'nag9', 'fgm', 'sst', 'sst_bound', 'fgm_bound')}
    for N in Ns:
        res['gd'].append(pep(cyclic_gd([1.0] * N)))
        res['smm'].append(pep(smm(N)))
        res['nag585'].append(pep(nag_const(N, BETA)))
        res['nag9'].append(pep(nag_const(N, 0.9)))
        Hf, Af = similar_triangles(N, [True] * N); res['fgm'].append(pep(Hf)); res['fgm_bound'].append(1 / (2 * Af))
        Hs, As = similar_triangles(N, [t % 2 == 0 for t in range(N)]); res['sst'].append(pep(Hs)); res['sst_bound'].append(1 / (2 * As))
        print(N, ' '.join(f"{k}={res[k][-1]:.5f}" for k in res), flush=True)
    res['N'] = Ns
    json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'ch5d_pep.json'), 'w'))
