import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
# Gain weight theorem (Chapter 7): verification by exact periodic Lyapunov equations.
# A linear periodic first order method on the scalar quadratic f = lambda e^2/2 with noisy
# gradient g = lambda * (evaluation point error) + xi.  At residue r:
#   y_t = gamma_r . z_t,  e_{t+1} = delta_r . z_t - h_r (a * y_t + alpha xi_t),  z_t = (e_t, ..., e_{t-K+1}).
# Units: alpha = sigma = 1, so a = lambda.
import numpy as np
from scipy.linalg import solve_discrete_lyapunov
import json, itertools

K = 3   # memory length used for all schedules (padding with zeros)

def pad(v):
    v = list(v) + [0.0] * (K - len(v)); return np.array(v, float)

def step_matrix(delta, gamma, h, a):
    S = np.zeros((K, K)); S[0, :] = delta - a * h * gamma
    for i in range(1, K): S[i, i - 1] = 1.0
    return S

def u():
    v = np.zeros(K); v[0] = 1.0; return v

def exact_variance(sched, a):
    """phase averaged stationary variance of e_t, and per step spectral radius"""
    p = len(sched)
    Ss = [step_matrix(d, g, h, a) for (d, g, h) in sched]
    A = np.eye(K); W = np.zeros((K, K))
    for r in range(p):
        A = Ss[r] @ A
        W = Ss[r] @ W @ Ss[r].T + sched[r][2] ** 2 * np.outer(u(), u())
    rad = max(abs(np.linalg.eigvals(A)))
    if rad >= 1: return np.inf, rad ** (1 / p)
    Sig = solve_discrete_lyapunov(A, W)
    vs = []
    for r in range(p):
        vs.append(Sig[0, 0]); Sig = Ss[r] @ Sig @ Ss[r].T + sched[r][2] ** 2 * np.outer(u(), u())
    return np.mean(vs), rad ** (1 / p)

def gains_by_impulse(sched, periods=4000):
    """c_r = minus the eventual displacement caused by a unit gradient impulse at residue r, with a = 0"""
    p = len(sched); cs = []
    for r in range(p):
        z = np.zeros(K)
        for t in range(p * periods):
            d, g, h = sched[t % p]
            new = d @ z - (h if t == r else 0.0)
            z = np.concatenate([[new], z[:-1]])
        cs.append(-z[0])
    return np.array(cs)

def gains_by_eigvec(sched):
    p = len(sched); Ss = [step_matrix(d, g, h, 0.0) for (d, g, h) in sched]
    A = np.eye(K)
    for S in Ss: A = S @ A
    w, V = np.linalg.eig(A.T); i = np.argmin(abs(w - 1)); l = np.real(V[:, i]); l = l / l.sum()
    others = sorted(abs(np.delete(w, i)))
    cs = []
    for r in range(p):
        P = np.eye(K)
        for s in range(r + 1, p): P = Ss[s] @ P
        cs.append(sched[r][2] * l @ P @ u())
    return np.array(cs), (others[-1] if others else 0.0)

D = (pad([1]), pad([1]), 1.0)                       # plain gradient step
def NAG(b): return (pad([1 + b, -b]), pad([1 + b, -b]), 1.0)
def HB(b): return (pad([1 + b, -b]), pad([1]), 1.0)
def GDh(h): return (pad([1]), pad([1]), h)

b0 = np.log(1.5) / np.log(2)
schedules = {
    'GD': [D],
    'HB, every, b0': [HB(b0)],
    'NAG, every, b0': [NAG(b0)],
    'NAG, every, 0.9': [NAG(0.9)],
    'SMM (NAG odd), b0': [D, NAG(b0)],
    'SMM (NAG odd), 0.9': [D, NAG(0.9)],
    'HB odd, b0': [D, HB(b0)],
    'cyclic GD (1,3)': [GDh(1.0), GDh(3.0)],
    'p3 NAG at {2}, b0': [D, D, NAG(b0)],
    'p3 NAG at {1,2}, b0': [D, NAG(b0), NAG(b0)],
    'p3 HB at {1,2}, 0.9': [D, HB(0.9), HB(0.9)],
    'p4 NAG at {3}, b0': [D, D, D, NAG(b0)],
    'p4 NAG at {2,3}, 0.9': [D, D, NAG(0.9), NAG(0.9)],
    'p4 NAG at {1,2,3}, b0': [D, NAG(b0), NAG(b0), NAG(b0)],
    'p4 NAG at {1,2,3}, 0.9': [D, NAG(0.9), NAG(0.9), NAG(0.9)],
}
# a generic schedule with unrelated memory and extrapolation coefficients and three term memory
rng = np.random.default_rng(7)
gen = [(pad([1.2, -0.3, 0.1]), pad([1.5, -0.5]), 0.8), (pad([0.9, 0.1]), pad([0.7, 0.3]), 1.3),
       (pad([1.4, -0.6, 0.2]), pad([1.0]), 0.6)]
schedules['generic p3 (three term memory)'] = gen

rows = []
for name, sch in schedules.items():
    p = len(sch)
    c_imp = gains_by_impulse(sch); c_eig, fast = gains_by_eigvec(sch)
    M = c_imp.sum(); S2 = (c_imp ** 2).sum()
    F = p * S2 / M ** 2
    res = {'name': name, 'p': p, 'c': [round(x, 6) for x in c_imp], 'max|c_imp-c_eig|': float(abs(c_imp - c_eig).max()),
           'fast_radius': float(fast), 'M': M, 'speed M/p': M / p, 'floor factor': F}
    for a in [1e-2, 1e-3, 1e-4]:
        v, r = exact_variance(sch, a)
        lead = S2 / (2 * M * a)
        res[f'var/lead a={a}'] = v / lead
        res[f'(1-r)/(Ma/p) a={a}'] = (1 - r) / (M * a / p)
        # matched speed floor relative to GD with step M/p: GD variance = (M/p)/(a(2 - a M/p)) exactly
        hg = M / p; vgd = hg / (a * (2 - a * hg))
        res[f'exact matched ratio a={a}'] = v / vgd
    rows.append(res)
    print(f"{name:32s} p={p} c={np.round(c_imp,4)} M/p={M/p:.4f} F={F:.4f} fast={fast:.3f} "
          f"var/lead(1e-3)={res['var/lead a=0.001']:.4f} var/lead(1e-4)={res['var/lead a=0.0001']:.5f} "
          f"speed(1e-4)={res['(1-r)/(Ma/p) a=0.0001']:.5f} matched(1e-3)={res['exact matched ratio a=0.001']:.4f} "
          f"eigcheck={res['max|c_imp-c_eig|']:.1e}")

# closed form for momentum on a subset: c_r = sum_{j=0}^{n_r} b^j, n_r = run of momentum residues after r
def subset_gains(mask, b):
    p = len(mask); cs = []
    for r in range(p):
        if all(mask): cs.append(1 / (1 - b)); continue
        n = 0; s = (r + 1) % p
        while mask[s]: n += 1; s = (s + 1) % p
        cs.append((1 - b ** (n + 1)) / (1 - b))
    return np.array(cs)

print("\nsubset schedules: closed form gains versus impulse gains, and best proper subset per p")
best = {}
allpts = []
for b in [0.5, b0, 0.9]:
    for p in range(1, 7):
        for mask in itertools.product([0, 1], repeat=p):
            sch = [NAG(b) if m else D for m in mask]
            c = gains_by_impulse(sch, periods=3000) if p <= 3 else None
            cf = subset_gains(mask, b)
            if c is not None: assert abs(c - cf).max() < 1e-8, (mask, c, cf)
            F = p * (cf ** 2).sum() / cf.sum() ** 2; spd = cf.sum() / p
            allpts.append((b, p, mask, spd, F))
            if 0 < sum(mask) < p:
                key = (round(b, 3), p)
                if key not in best or F < best[key][1]: best[key] = (mask, F, spd)
for k, v in sorted(best.items()): print(k, 'best proper subset', v[0], 'F=%.4f speed=%.4f' % (v[1], v[2]))
json.dump({'rows': rows, 'subset_points': [(b, p, list(m), s, F) for (b, p, m, s, F) in allpts]},
          open(BOOK + '/scripts/ch7_gainweight.json', 'w'), indent=1, default=float)

# rows for the book table: exact matched speed floor ratio at a = 1e-2 and 1e-3 against the predicted F
print("\nTABLE")
for res in rows:
    print(f"{res['name']:32s} c={res['c']} M/p={res['speed M/p']:.4f} F={res['floor factor']:.4f} "
          f"exact(1e-2)={res['exact matched ratio a=0.01']:.4f} exact(1e-3)={res['exact matched ratio a=0.001']:.4f} "
          f"speed ratio(1e-3)={res['(1-r)/(Ma/p) a=0.001']:.4f} var/lead(1e-3)={res['var/lead a=0.001']:.4f}")
# single injection per period: F = 1 + (p-1) b^2/(p+b)^2
for p in range(2, 7):
    print(p, 'single injection F formula', 1 + (p - 1) * b0 ** 2 / (p + b0) ** 2)
print("\nTABLE2 exact matched ratio at a = 1e-2, 1e-3, 1e-4")
for res in rows:
    print(f"{res['name']:32s} F={res['floor factor']:.4f} " + " ".join(f"{res[f'exact matched ratio a={a}']:.4f}" for a in [0.01, 0.001, 0.0001]))
