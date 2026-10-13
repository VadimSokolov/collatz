"""Chapter 8: ladder thresholds (P2), Diophantine insufficiency numbers, and the
empirical irrationality profile of log_2 3.  Figures ch8_thresholds.pdf and
ch8_irrprofile.pdf.  High precision via mpmath; continued fraction confirmed at
two precisions."""
import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
import sys; sys.path.insert(0, BOOK + '/scripts'); from bookstyle import *
import mpmath as mp, numpy as np, math
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def cf_of(dps, n):
    mp.mp.dps = dps
    th = mp.log(3)/mp.log(2); y = th; a = []
    for _ in range(n):
        ai = int(mp.floor(y)); a.append(ai); y = 1/(y-ai)
    return a
A1 = cf_of(700, 260); A2 = cf_of(500, 260)
NCF = next(i for i in range(260) if A1[i] != A2[i]) if A1 != A2 else 260
NCF = min(NCF - 5, 200)
A = A1[:NCF]
print('partial quotients confirmed at two precisions:', NCF)
mp.mp.dps = 700
TH = mp.log(3)/mp.log(2); L2 = mp.log(2)
p = [0, 1]; q = [1, 0]          # p_{-2}, p_{-1}; q_{-2}, q_{-1}
for ai in A:
    p.append(ai*p[-1]+p[-2]); q.append(ai*q[-1]+q[-2])
P = lambda k: p[k+2]; Q = lambda k: q[k+2]
eps = lambda k: abs(Q(k)*TH-P(k))
delta = lambda r: mp.ceil(r*TH)-r*TH
F = lambda r: r/(3*delta(r)*L2)

# ---------------- upper ladder rungs (delta records) ----------------
rungs = []
for k in range(0, NCF-2, 2):
    for j in range(1, A[k+1]+1):
        rungs.append((Q(k-1)+j*Q(k), k, j))
rungs.sort()
LIM = 10**18
rungs = [t for t in rungs if t[0] <= LIM]
Fr = [(u, k, j, F(u)) for (u, k, j) in rungs]
# check F increasing along rungs up to 1e18
inc = all(Fr[i+1][3] > Fr[i][3] for i in range(len(Fr)-1))
print('rungs <= 1e18:', len(Fr), ' F increasing along all of them:', inc)
bad = [(Fr[i][0], Fr[i+1][0]) for i in range(len(Fr)-1) if not Fr[i+1][3] > Fr[i][3]]
print('non increasing pairs (first 10):', bad[:10])
print('number of rungs <= 72057431991:', sum(1 for t in Fr if t[0] <= 72057431991))

# ---------------- P2: thresholds ----------------
# r_*(X) = least rung u with F(u) > X.  It jumps at X = running maxima of F.
print('\nThresholds: r_*(X) for log2 X in [60, 90]')
runmax = mp.mpf(0); rec = []
for (u, k, j, f) in Fr:
    if f > runmax:
        rec.append((u, k, j, f, runmax)); runmax = f
# For X in [F(prev record), F(u)), r_* = u.
rows = []
for (u, k, j, f, prev) in rec:
    lo = float(mp.log(prev, 2)) if prev > 0 else -1
    hi = float(mp.log(f, 2))
    if hi >= 60 and lo <= 90:
        rows.append((u, int(mp.ceil(u*TH)), k, j, lo, hi))
        print(f'  r_* = {u:>22d}  s = {int(mp.ceil(u*TH)):>22d}  rung k={k} j={j}   for {lo:7.3f} <= log2 X < {hi:7.3f}')
HERCHER_X = 1536*2**60; HERCHER_R = 1.375e11
print('log2 of Hercher verification range 1536*2^60 =', math.log2(HERCHER_X))

# ---------------- Diophantine insufficiency ----------------
print('\nInsufficiency proposition numbers')
# min over rungs of r*delta(r) gives the best constant c in delta(r) >= c/r (one sided mu = 2)
for lim in [72057431991, 10**13, LIM]:
    cands = [(u*delta(u), u) for (u, k, j, f) in Fr if u <= lim]
    cmin, uarg = min(cands)
    print(f'  min r delta(r) over rungs <= {lim:.3g}: {mp.nstr(cmin, 6)} at r = {uarg}')
cobs = min((u*delta(u), u) for (u, k, j, f) in Fr if u <= 10**13)[0]
for e in [40, 60, 68, 72, 80]:
    X = mp.mpf(2)**e
    rmu2 = mp.sqrt(3*cobs*X*L2)
    epsX = 1/(3*X*L2)
    rstar = next(u for (u, k, j, f) in Fr if f > X)
    print(f'  X=2^{e}: mu=2 measure (c={mp.nstr(cobs,4)}) excludes r < {mp.nstr(rmu2, 5)};  exact r_* = {rstar};  ratio {mp.nstr(rstar/rmu2, 4)}')
# convergents from above that violate the exclusion at X = 2^68: first few
X = mp.mpf(2)**68; epsX = 1/(3*X*L2)
viol = [Q(k) for k in range(1, 60, 2) if delta(Q(k)) < Q(k)*epsX][:5]
print('  first odd convergent denominators with delta < r eps at 2^68:', viol)
# Best conceivable constants: if |theta - s/r| >= c r^-mu for all r, then c <= c_mu :=
# min_r r^(mu-1) delta(r); the minimum over all r is attained at a delta record (rung).
print('  best conceivable constants c_mu and the r excluded at X (r^mu <= 3 c_mu X ln 2):')
for mu in [2, 2.5, 3, 4, 5]:
    cm_, um_ = min(((mp.mpf(u)**(mu-1))*delta(u), u) for (u, k, j, f) in Fr)
    ex = {e: (3*cm_*mp.mpf(2)**e*L2)**(1/mp.mpf(mu)) for e in (60, 68, 72)}
    print(f'   mu={mu}: c_mu <= {mp.nstr(cm_,5)} (attained at r={um_});  excluded r up to: '
          + ', '.join(f'2^{e}: {mp.nstr(v,4)}' for e, v in ex.items()))
# check q_k q_{k+1} >= 3 X ln 2 criterion for odd k at 2^68
crit = [Q(k) for k in range(1, 40, 2) if Q(k)*Q(k+1) >= 3*X*L2][:3]
print('  odd k with q_k q_(k+1) >= 3 X ln 2 at 2^68, first three q_k:', crit)

# ---------------- irrationality profile ----------------
print('\nIrrationality profile nu_k = -log|theta - p_k/q_k| / log q_k')
K = list(range(2, NCF-1))
muk = []
for k in K:
    err = abs(TH - mp.mpf(P(k))/Q(k))
    muk.append(float(-mp.log(err)/mp.log(Q(k))))
for k, m in zip(K[:30], muk[:30]):
    print(f'  k={k:3d} a_(k+1)={A[k+1]:4d} q_k={mp.nstr(mp.mpf(Q(k)),4):>10s}  nu_k={m:.4f}')
print('  max nu_k for k>=10:', max(muk[8:]), ' at k =', K[8+int(np.argmax(muk[8:]))])
print('  mean nu_k for k>=50:', np.mean(muk[48:]), ' last k:', K[-1], 'mu:', muk[-1])
print('  largest partial quotients among a_1..a_%d:' % (NCF-1), sorted(((a, i) for i, a in enumerate(A) if i > 0), reverse=True)[:8])
# Khinchin's constant check (geometric mean of partial quotients)
gm = math.exp(np.mean([math.log(a) for a in A[1:]]))
print('  geometric mean of a_1..a_n:', gm)

# ---------------- figure: thresholds step plot ----------------
fig, ax = plt.subplots(figsize=(WIDTH, 2.9))
xs = []; ys = []
for (u, s, k, j, lo, hi) in rows:
    lo2 = max(lo, 60); hi2 = min(hi, 90)
    xs += [lo2, hi2]; ys += [u, u]
ax.fill_between(xs, ys, 1e8, color=PAL['sky'], alpha=0.18, lw=0)
ax.plot(xs, ys, color=PAL['teal'], lw=1.6, label=r'sharp bound $r_*(X)$')
for (u, s, k, j, lo, hi) in rows:
    if 60 <= lo <= 90:
        ax.plot(lo, u, 'o', ms=3.2, color=PAL['teal'], mfc='white', mew=1.0)
ax.axvline(68, color=PAL['coral'], lw=1.1, ls='--')
ax.text(66.0, 3e13, 'Ba\u0159ina\n$2^{68}$', color=PAL['coral'], fontsize=7, va='top', ha='right')
ax.axvline(71, color=PAL['coral'], lw=1.1, ls=':')
ax.text(71.4, 3e13, 'Ba\u0159ina\n$2^{71}$', color=PAL['coral'], fontsize=7, va='top')
hx = math.log2(HERCHER_X)
ax.plot(hx, HERCHER_R, marker='*', ms=10, color=PAL['gold'], mec=PAL['ink'], mew=0.5, ls='none',
        label=r'Hercher, $r\geq1.375\times10^{11}$ at $1536\cdot2^{60}$')
ax.set_yscale('log'); ax.set_xlim(60, 90); ax.set_ylim(1e8, 3e15)
ax.set_xlabel(r'verified range $\log_2 X$')
ax.set_ylabel(r'least possible $r$')
ax.legend(loc='lower right', fontsize=7)
fig.savefig(BOOK + '/figures/ch8_thresholds.pdf'); plt.close(fig)

# ---------------- figure: irrationality profile ----------------
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(WIDTH, 3.6), sharex=True,
                               gridspec_kw=dict(height_ratios=[2.2, 1], hspace=0.08))
a_next = np.array([A[k+1] for k in K])
sc = ax1.scatter(K, muk, c=a_next, cmap=CMAP_HOT, norm=LogNorm(1, a_next.max()), s=12,
                 edgecolor=PAL['ink'], linewidth=0.25, zorder=3)
ax1.plot(K, muk, color=PAL['grey'], lw=0.5, zorder=2)
ax1.axhline(2, color=PAL['teal'], lw=1.0, ls='--')
ax1.text(K[-1], 2.02, r'$\nu=2$', color=PAL['teal'], ha='right', va='bottom', fontsize=7)
ax1.set_ylabel(r'$\nu_k$')
cb = fig.colorbar(sc, ax=[ax1, ax2], pad=0.02, aspect=30)
cb.set_label(r'next partial quotient $a_{k+1}$')
ax2.bar(K, a_next, color=[plt.get_cmap(CMAP_HOT)(LogNorm(1, a_next.max())(v)) for v in a_next],
        width=0.85, edgecolor=PAL['ink'], linewidth=0.2)
ax2.set_yscale('log'); ax2.set_ylabel(r'$a_{k+1}$'); ax2.set_xlabel(r'convergent index $k$')
ax2.grid(axis='x', visible=False)
fig.savefig(BOOK + '/figures/ch8_irrprofile.pdf'); plt.close(fig)
