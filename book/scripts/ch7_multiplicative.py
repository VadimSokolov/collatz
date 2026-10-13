import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
# Chapter 7, multiplicative (random curvature) noise.  Scalar eigendirection; the gradient at evaluation point y is
#   g = (lambda + eta_t) e_y + xi_t,  E eta = E xi = E eta xi = 0,  Var(alpha eta_t) = omega^2,  Var xi_t = sigma^2,
# fresh (eta_t, xi_t) at each step, independent of the past.  Units alpha = sigma = 1 so a = alpha lambda.
# Least squares scaling omega^2 = theta a^2 (theta = Var(x^2)/(E x^2)^2 for a scalar regressor x; 2 if Gaussian).
import numpy as np, sympy as sp
b0 = np.log(1.5) / np.log(2)

# ---------- 1. exact SMM second moments (symbolic) ----------
a, b, w2, th, eps = sp.symbols('a beta omega2 theta epsilon', positive=True)
c1 = (1 - a) ** 2 + w2
c2 = (1 - (1 + b) * a) ** 2 + (1 + b) ** 2 * w2
me = (c1 * (1 + b) ** 2 + 1) / (1 - c1 * c2)
mo = c1 * me + 1
msmm = (me + mo) / 2
mgd = 1 / (1 - c1)
V = (2 + 2 * b + b ** 2) / (2 + b)
print("[1] least squares scaling: lim a*m_SMM =", sp.simplify(sp.limit(a * msmm.subs(w2, th * a ** 2), a, 0)),
      "; lim ratio SMM/GD =", sp.simplify(sp.limit((msmm / mgd).subs(w2, th * a ** 2), a, 0)))
lg = sp.limit(a * mgd.subs(w2, eps * a), a, 0); ls_ = sp.limit(a * msmm.subs(w2, eps * a), a, 0)
print("[1] diffusive scaling omega^2 = eps a: a*m_GD ->", sp.simplify(lg), "; a*m_SMM ->", sp.simplify(ls_))
pred = ((1 + b) ** 2 + 1) / (2 * (2 + b) - eps * ((1 + b) ** 2 + 1))
print("    matches Sum c^2/(2M - eps Sum c^2):", sp.simplify(ls_ - pred) == 0)

# ---------- 2. general periodic linear method: exact second moment operator ----------
K = 3
def pad(v): v = list(v) + [0.0] * (K - len(v)); return np.array(v, float)
def comp(delta):
    C = np.zeros((K, K)); C[0, :] = delta
    for i in range(1, K): C[i, i - 1] = 1.0
    return C
u = np.zeros(K); u[0] = 1.0
D = (pad([1]), pad([1]), 1.0)
def NAG(bb): return (pad([1 + bb, -bb]), pad([1 + bb, -bb]), 1.0)
def HB(bb): return (pad([1 + bb, -bb]), pad([1]), 1.0)
def step_op(sch_r, aa, om2):
    d, g, h = sch_r; S = comp(d) - aa * h * np.outer(u, g)
    G = np.outer(u, g)
    I = np.eye(K)
    # vec(S X S^T) = (S kron S) vec X ; vec(G X G^T) = (G kron G) vec X
    return np.kron(S, S) + om2 * h * h * np.kron(G, G), h * h * np.outer(u, u).reshape(-1)
def exact_m2(sch, aa, om2):
    """phase averaged stationary E e^2 (units alpha^2 sigma^2), or inf if not mean square stable"""
    n = K * K; Tm = np.eye(n); c = np.zeros(n); ops = []
    for r in range(len(sch)):
        Op, src = step_op(sch[r], aa, om2); ops.append((Op, src)); Tm = Op @ Tm; c = Op @ c + src
    if max(abs(np.linalg.eigvals(Tm))) >= 1 - 1e-13: return np.inf
    X = np.linalg.solve(np.eye(n) - Tm, c); vals = []
    for Op, src in ops:
        vals.append(X.reshape(K, K)[0, 0]); X = Op @ X + src
    return float(np.mean(vals))
def num(expr, aa, bb, om2):
    return float(expr.subs({a: aa, b: bb, w2: om2}))
for (aa, om2) in [(0.05, 0.01), (0.3, 0.05), (0.8, 0.3)]:
    print(f"[2] check a={aa} omega2={om2}: SMM closed {num(msmm, aa, b0, om2):.6f} operator {exact_m2([D, NAG(b0)], aa, om2):.6f};"
          f" GD closed {num(mgd, aa, b0, om2):.6f} operator {exact_m2([D], aa, om2):.6f}")

# ---------- 3. Monte Carlo check (Rademacher curvature noise, Gaussian additive noise) ----------
rng = np.random.default_rng(2026)
def mc(kind, aa, bb, om2, R=4000, T=6000, burn=2000):
    x = np.zeros(R); xp = np.zeros(R); acc = 0.0; n = 0
    for t in range(T):
        at = aa + np.sqrt(om2) * rng.choice([-1.0, 1.0], R); xi = rng.standard_normal(R)
        mom = (kind == 'nag') or (kind == 'smm' and t % 2 == 1)
        y = x + bb * (x - xp) if mom else x
        xp, x = x, y - (at * y + xi)
        if t >= burn: acc += np.mean(x * x); n += 1
    return acc / n
for (aa, th_) in [(0.1, 1.0), (0.3, 0.5)]:
    om2 = th_ * aa * aa
    print(f"[3] a={aa} theta={th_}: SMM exact {exact_m2([D, NAG(b0)], aa, om2):.4f} MC {mc('smm', aa, b0, om2):.4f};"
          f" NAG exact {exact_m2([NAG(b0)], aa, om2):.4f} MC {mc('nag', aa, b0, om2):.4f}")

# ---------- 4. diffusive scaling: periodic gain prediction versus exact ----------
def gains(sch):
    # impulse response at a = 0
    p = len(sch); cs = []
    for r in range(p):
        z = np.zeros(K)
        for t in range(p * 3000):
            d, g, h = sch[t % p]; new = d @ z - (h if t == r else 0.0); z = np.concatenate([[new], z[:-1]])
        cs.append(-z[0])
    return np.array(cs)
print("[4] diffusive scaling omega^2 = eps a: exact a*m versus Sum c^2/(2M - eps Sum c^2), and threshold 2M/Sum c^2")
for name, sch in [('GD', [D]), ('SMM b0', [D, NAG(b0)]), ('NAG b0', [NAG(b0)]), ('HB b0', [HB(b0)]),
                  ('SMM 0.9', [D, NAG(0.9)]), ('NAG 0.9', [NAG(0.9)]), ('p3 HB {1,2} 0.9', [D, HB(0.9), HB(0.9)])]:
    c = gains(sch); M = c.sum(); S2 = (c ** 2).sum(); thr = 2 * M / S2
    out = []
    for e_ in [0.25 * thr, 0.5 * thr, 0.8 * thr]:
        aa = 1e-4; ex = aa * exact_m2(sch, aa, e_ * aa); pr = S2 / (2 * M - e_ * S2); out.append(f"{ex / pr:.4f}")
    print(f"    {name:16s} c={np.round(c, 4)} threshold eps* = {thr:.4f}; exact/pred at eps = (0.25, 0.5, 0.8) eps*: {out}")
print("    SMM threshold 2/V(beta) at b0:", 2 / float(V.subs(b, b0)), " NAG threshold 2(1-beta):", 2 * (1 - b0))

# ---------- 5. least squares scaling: stability limits and the uniform bound ----------
def stab_limit(sch, th_, amax=2.5):
    lo, hi = 1e-6, amax
    A = np.linspace(1e-3, amax, 2500); ok = [np.isfinite(exact_m2(sch, x, th_ * x * x)) for x in A]
    if all(ok): return amax
    i = ok.index(False); lo, hi = A[i - 1], A[i]
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if np.isfinite(exact_m2(sch, mid, th_ * mid * mid)): lo = mid
        else: hi = mid
    return lo
print("[5] mean square stability limits in a (largest stable a), least squares scaling")
TH = [0.0, 0.25, 0.5, 1.0, 2.0]
lims = {}
for name, sch in [('GD', [D]), ('SMM b0', [D, NAG(b0)]), ('SMM 0.9', [D, NAG(0.9)]), ('NAG b0', [NAG(b0)]), ('NAG 0.9', [NAG(0.9)])]:
    lims[name] = [stab_limit(sch, t) for t in TH]
    print(f"    {name:8s}", " ".join(f"{x:.4f}" for x in lims[name]))
for bb in [b0, 0.9]:
    tmax = (-bb ** 2 + np.sqrt(bb ** 4 + 4 * (1 + bb) ** 2)) / (2 * (1 + bb) ** 2)
    print(f"    SMM stable at a = 1 iff theta < {tmax:.4f} (beta = {bb:.3f})")
print("[5] sup over a in (0, 1] (both stable) of m/m_GD; SMM compared with V(beta)")
sup = {}
Agrid = np.concatenate([np.logspace(-3, -1, 60), np.linspace(0.1, 1.0, 400)])
for bb in [b0, 0.9]:
    Vb = float(V.subs(b, bb))
    for t in TH:
        rs, rn = [], []
        for x in Agrid:
            g_ = exact_m2([D], x, t * x * x)
            rs.append(exact_m2([D, NAG(bb)], x, t * x * x) / g_ if np.isfinite(g_) else np.nan)
            rn.append(exact_m2([NAG(bb)], x, t * x * x) / g_ if np.isfinite(g_) else np.nan)
        rs, rn = np.array(rs), np.array(rn); sup[(round(bb, 3), t)] = (rs, rn)
        fs = np.isfinite(rs); fn = np.isfinite(rn)
        msg = f"    beta={bb:.3f} theta={t}: SMM sup ratio {np.nanmax(np.where(fs, rs, np.nan)):.4f} (V={Vb:.4f})"
        if (~fs & np.isfinite([exact_m2([D], x, t * x * x) for x in Agrid])).any(): msg += " SMM unstable somewhere in (0,1]"
        msg += f"; NAG sup ratio {np.nanmax(np.where(fn, rn, np.nan)):.4f} (slow {1/(1-bb):.3f})"
        if (~fn).any(): msg += f" NAG unstable for a >= {Agrid[~fn][0]:.3f}"
        print(msg)

# ---------- 6. figure ----------
import sys; sys.path.insert(0, BOOK + '/scripts'); from bookstyle import *
fig, axs = plt.subplots(2, 1, figsize=(WIDTH, 5.0), gridspec_kw=dict(hspace=0.45))
ax = axs[0]
cols = {0.0: PAL['ink'], 0.5: PAL['teal'], 1.0: PAL['violet'], 2.0: PAL['coral']}
# Grayscale: teal and violet print alike, so the marker shape tells the noise level; solid and dashed still tell SMM from NAG.
mkth = {0.0: 'o', 0.5: 's', 1.0: '^', 2.0: 'D'}
def mk_at(y, ph, n=5, lo=0.1, hi=4.3):
    ok = np.where((Agrid >= lo) & np.isfinite(y) & (y < hi))[0]
    return list(ok[np.minimum(((np.arange(n) + ph) / n * len(ok)).astype(int), len(ok) - 1)])
hth = []
for i, (t, c) in enumerate(cols.items()):
    rs, rn = sup[(round(b0, 3), t)]
    ax.plot(Agrid, rs, color=c, lw=2.0, marker=mkth[t], ms=4.4, mec='white', mew=0.5, markevery=mk_at(rs, 0.15 + 0.12 * i))
    ax.plot(Agrid, rn, color=c, lw=1.1, ls='--', marker=mkth[t], ms=3.8, mec='white', mew=0.4, markevery=mk_at(rn, 0.6 + 0.12 * i))
    hth.append(plt.Line2D([], [], color=c, lw=2.0, marker=mkth[t], ms=4.4, mec='white', mew=0.5, label=rf'$\vartheta={t:g}$'))
Vb0 = float(V.subs(b, b0))
ax.axhline(Vb0, color=PAL['gold'], lw=1.0, ls=':'); ax.text(1.15e-3, Vb0 * 1.05, r'$V(\beta)=1.359$', fontsize=7.5, color=PAL['ink'])
ax.axhline(1 / (1 - b0), color=PAL['grey'], lw=1.0, ls=':'); ax.text(1.15e-3, 2.409 * 1.04, r'$1/(1-\beta)=2.409$', fontsize=7.5, color=PAL['ink'])
ax.set_xscale('log'); ax.set_xlim(1e-3, 1); ax.set_ylim(0.8, 4.6)
ax.set_xlabel(r'normalized step $a=\alpha\lambda$'); ax.set_ylabel('ratio to GD')
ax.legend(handles=hth, loc='upper left', bbox_to_anchor=(0.0, 0.93), fontsize=7, ncol=4, columnspacing=0.8, handlelength=1.8,
          title=r'solid SMM, dashed NAG, $\beta=0.585$', title_fontsize=7)
ax.set_title('(a)', loc='left', fontsize=9)
ax = axs[1]
T2 = np.linspace(0, 2.0, 81)
for (name, sch, lab, c, lsty) in [('GD', [D], 'GD', PAL['green'], (0, (3.2, 1.2, 0.8, 1.2))),   # dash-dot: GD and SMM beta=0.9 print alike in grayscale
                                  ('SMM b0', [D, NAG(b0)], r'SMM, $\beta=0.585$', PAL['gold'], '-'),
                                  ('SMM 0.9', [D, NAG(0.9)], r'SMM, $\beta=0.9$', PAL['coral'], '-'),
                                  ('NAG b0', [NAG(b0)], r'NAG, $\beta=0.585$', PAL['violet'], '--'),
                                  ('NAG 0.9', [NAG(0.9)], r'NAG, $\beta=0.9$', PAL['sky'], '--')]:
    ax.plot(T2, [stab_limit(sch, t) for t in T2], color=c, lw=1.8, ls=lsty, label=lab)
ax.axhline(1.0, color=PAL['grey'], lw=0.8, ls=':'); ax.text(1.62, 1.04, r'$a=1$ ($\alpha=1/L$)', fontsize=7.5, color=PAL['ink'])
ax.set_xlabel(r'relative curvature noise $\vartheta$'); ax.set_ylabel(r'largest mean square stable $a$')
ax.set_ylim(0, 2.1); ax.set_xlim(0, 2); ax.legend(loc='upper right', fontsize=7, ncol=2, handlelength=2.8)
ax.set_title('(b)', loc='left', fontsize=9)
fig.savefig(BOOK + '/figures/ch7_multiplicative.pdf')
