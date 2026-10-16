"""Figures for the chapter on three point inequalities and similar triangles.
ch5d_geometry.pdf : (a) the three point identity as the law of cosines; (b) similar triangles.
ch5d_pep.pdf      : (a) exact worst cases on smooth convex functions; (b) the threshold in beta.
ch5d_collatz.pdf  : (a) Bregman decomposition along one orbit; (b) dissipation values and frequencies."""
import sys, os, json, math, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bookstyle import *
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); FIG = os.path.join(HERE, '..', 'figures'); DATA = os.path.join(HERE, 'data')
BETA = math.log(1.5) / math.log(2); PHI1 = (math.sqrt(5) - 1) / 2

# ---------------------------------------------------------------- geometry
fig, axs = plt.subplots(1, 2, figsize=(WIDTH, 2.25), gridspec_kw=dict(wspace=0.08))
ax = axs[0]; ax.set_aspect('equal'); ax.axis('off')
u = np.array([0.15, 0.0]); z = np.array([2.75, 1.3]); zp = np.array([1.75, 0.15])
for a, b, c in ((u, z, PAL['ink']), (u, zp, PAL['ink']), (z, zp, PAL['coral'])):
    ax.plot([a[0], b[0]], [a[1], b[1]], color=c, lw=1.2)
for p, t, off in ((u, '$u$', (-0.05, -0.28)), (z, '$z$', (0.06, 0.05)), (zp, '$z^{+}$', (0.08, -0.24))):
    ax.plot(*p, 'o', color=PAL['ink'], ms=3.5); ax.text(p[0] + off[0], p[1] + off[1], t, fontsize=8.5)
g = (zp - z); ax.annotate('', xy=zp + 0.0, xytext=z, arrowprops=dict(arrowstyle='-|>', color=PAL['coral'], lw=1.1))
ax.text(2.35, 0.55, r'$-\eta g$', color=PAL['coral'], fontsize=8)
ax.text(-0.1, 2.15, '(a) three point identity', fontsize=7.5)
ax.text(-0.1, 1.88, r'$2\langle z-z^{+},\,u-z^{+}\rangle$', fontsize=7, color=PAL['ink'])
ax.text(-0.1, 1.6, r'$=|u-z^{+}|^2+|z-z^{+}|^2-|u-z|^2$', fontsize=7, color=PAL['ink'])
ax.set_xlim(-0.3, 3.1); ax.set_ylim(-0.45, 2.35)
ax = axs[1]; ax.set_aspect('equal'); ax.axis('off')
x = np.array([0.0, 0.0]); zk = np.array([2.7, 1.7]); zk1 = np.array([2.2, 0.2]); tau = 0.42
y = x + tau * (zk - x); x1 = x + tau * (zk1 - x)
ax.fill([x[0], zk[0], zk1[0]], [x[1], zk[1], zk1[1]], color=PAL['sky'], alpha=0.14, lw=0)
ax.fill([x[0], y[0], x1[0]], [x[1], y[1], x1[1]], color=PAL['gold'], alpha=0.35, lw=0)
for a, b, c, w in ((x, zk, PAL['ink'], 0.8), (x, zk1, PAL['ink'], 0.8), (zk, zk1, PAL['coral'], 1.2), (y, x1, PAL['teal'], 1.2)):
    ax.plot([a[0], b[0]], [a[1], b[1]], color=c, lw=w)
lab = ((x, '$x_k$', (-0.3, -0.15)), (zk, '$z_k$', (0.06, 0.04)), (zk1, '$z_{k+1}$', (0.06, -0.18)),
       (y, '$y_k$', (-0.3, 0.06)), (x1, '$x_{k+1}$', (-0.05, -0.27)))
for p, t, off in lab:
    ax.plot(*p, 'o', color=PAL['ink'], ms=3.2); ax.text(p[0] + off[0], p[1] + off[1], t, fontsize=8)
ax.text(2.55, 0.95, r'$-a\,g(y_k)$', color=PAL['coral'], fontsize=7.5)
ax.text(1.18, 0.5, r'$-g(y_k)$', color=PAL['teal'], fontsize=7.5)
ax.text(-0.25, 2.15, '(b) similar triangles', fontsize=7.5)
ax.text(-0.25, 1.85, r'ratio $\tau=a_{k+1}/A_{k+1}$', fontsize=7)
ax.set_xlim(-0.45, 3.2); ax.set_ylim(-0.4, 2.35)
fig.savefig(os.path.join(FIG, 'ch5d_geometry.pdf')); plt.close(fig)

# ---------------------------------------------------------------- PEP results
r = json.load(open(os.path.join(DATA, 'ch5d_pep.json'))); cr = json.load(open(os.path.join(DATA, 'ch5d_cross.json')))
N = np.array(r['N'] + [c[0] for c in cr])
cat = lambda k, idx: np.array(r[k] + [c[idx] for c in cr])
gd, sm, sst, fgm = cat('gd', 1), cat('smm', 2), cat('sst', 3), cat('fgm', 4)
fig, axs = plt.subplots(1, 2, figsize=(WIDTH, 2.35), gridspec_kw=dict(wspace=0.42))
ax = axs[0]; Nf = np.linspace(1, 40, 200)
ax.loglog(Nf, 1 / (4 * Nf + 2), color=PAL['grey'], lw=0.8, ls='--')
even = N % 2 == 0
ax.loglog(Nf, 1 / ((4 + 2 * BETA) * Nf + 2), color=PAL['teal'], lw=0.8, ls='--')
ax.loglog(N, gd, 'o', color=PAL['grey'], ms=2.6, label='gradient descent')
ax.loglog(N[even], sm[even], 's', color=PAL['teal'], ms=2.6, label='Switching MM')
Nn = np.array(r['N']); ax.loglog(Nn, r['nag585'], '^', color=PAL['violet'], ms=2.4, label=r'Nesterov, $\beta=0.585$')
ax.loglog(Nn, r['nag9'], 'v', color=PAL['plum'], ms=2.4, label=r'Nesterov, $\beta=0.9$')
ax.loglog(N, sst, 'D', color=PAL['gold'], ms=2.4, label='switched triangles')
ax.loglog(N, fgm, 'o', color=PAL['coral'], ms=2.4, label='similar triangles')
ax.set_xlabel('iterations $N$'); ax.set_ylabel(r'worst case of $(f(x_N)-f_*)/LR^2$')
ax.legend(fontsize=5.0, loc='lower left'); ax.set_title('(a)', loc='left', fontsize=8)
ax = axs[1]
th = {}
from ch5d_pep import pep, smm
betas = np.linspace(0.5, 0.98, 25)
for n, col in ((2, PAL['coral']), (4, PAL['violet']), (6, PAL['teal'])):
    vals = [pep(smm(n, b)) * ((4 + 2 * b) * n + 2) for b in betas]
    ax.plot(betas, vals, 'o-', color=col, ms=2.2, lw=0.9, label=f'$N={n}$')
b2 = json.load(open(os.path.join(DATA, 'ch5d_threshold.json')))['2']
ax.axvline(PHI1, color=PAL['gold'], lw=0.9, ls=':'); ax.text(PHI1 + 0.006, 1.86, r'$\varphi-1$', fontsize=6.0, color=PAL['gold'], rotation=90, va='top')
ax.axvline(b2, color=PAL['coral'], lw=0.9, ls='--'); ax.text(b2 + 0.008, 1.86, r'$\beta_2$', fontsize=6.5, color=PAL['coral'], va='top')
ax.axvline(BETA, color=PAL['ink'], lw=0.7); ax.text(BETA - 0.006, 1.86, 'Collatz', fontsize=6, ha='right', rotation=90, va='top')
ax.set_ylim(0.95, 1.9); ax.set_xlabel(r'$\beta$'); ax.set_ylabel('worst case / formula')
ax.legend(fontsize=5.5, loc='center right'); ax.set_title('(b)', loc='left', fontsize=8)
fig.savefig(os.path.join(FIG, 'ch5d_pep.pdf')); plt.close(fig)

# ---------------------------------------------------------------- Collatz decomposition
def S(v):
    w = 3 * v + 1
    while w % 2 == 0: w //= 2
    return w
def lg(v):
    b = v.bit_length(); s = max(0, b - 60); return s * math.log(2) + math.log(v >> s)
random.seed(11); v = random.getrandbits(3000) | 1
lin, dis, tot = [0.0], [0.0], [0.0]
while v > 1:
    w = S(v); lm = lg(w) - lg(v); M = math.exp(lm)
    lin.append(lin[-1] + M - 1); dis.append(dis[-1] - (M - 1 - lm)); tot.append(tot[-1] + lm); v = w
fig, axs = plt.subplots(1, 2, figsize=(WIDTH, 2.25), gridspec_kw=dict(wspace=0.42))
ax = axs[0]; t = np.arange(len(tot))
ax.plot(t, lin, color=PAL['violet'], lw=0.8, label=r'$\sum (M_i-1)$, martingale part')
ax.plot(t, dis, color=PAL['coral'], lw=0.9, label=r'$-\sum D(M_i)$, dissipation')
ax.plot(t, tot, color=PAL['ink'], lw=1.0, label=r'$\log y_n-\log y_0$')
ax.plot(t, -math.log(4 / 3) * t, color=PAL['grey'], lw=0.7, ls=':')
ax.set_xlabel('Syracuse steps $n$'); ax.legend(fontsize=5.2, loc='lower left'); ax.set_title('(a)', loc='left', fontsize=8)
ax = axs[1]
ks = np.arange(1, 9); Dk = [3 / 2 ** k - 1 - math.log(3 / 2 ** k) for k in ks]
ax.bar(ks, 2.0 ** -ks, width=0.6, color=PAL['sky'], alpha=0.6, label=r'model weight $2^{-k}$')
ax2 = ax.twinx(); ax2.plot(ks, Dk, 'o-', color=PAL['coral'], ms=3, lw=0.9, label=r'$D(3/2^k)$')
ax2.axhline(0.0291, color=PAL['ink'], lw=0.7, ls='--'); ax2.text(3.2, 0.034, r'$c_*=D(7/9)$', fontsize=6)
ax2.set_yscale('log'); ax2.set_ylabel('dissipation', fontsize=8); ax2.grid(False)
ax.set_xlabel('halving exponent $k$'); ax.set_ylabel('probability'); ax.set_title('(b)', loc='left', fontsize=8)
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels(); ax.legend(h1 + h2, l1 + l2, fontsize=5.5, loc='upper center', bbox_to_anchor=(0.62, 0.98)); ax.set_ylim(0, 0.75); ax2.set_ylim(0.015, 30)
fig.savefig(os.path.join(FIG, 'ch5d_collatz.pdf')); plt.close(fig)
print('ok')
