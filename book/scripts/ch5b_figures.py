import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
# Figures for the chapter on nonlinear difference equations and momentum MM.
import sys; sys.path.insert(0, BOOK + '/scripts')
from bookstyle import *
import numpy as np, json
from ch5b_jsr import mats, rho, jsr_exact, b0
FIG = BOOK + '/figures/'

# 1. Collatz as a switched system in log coordinates
fig, axs = plt.subplots(1, 2, figsize=(WIDTH, 2.2), gridspec_kw=dict(wspace=0.35))
ax = axs[0]
rng = np.random.default_rng(3)
for n0 in rng.integers(10**6, 10**7, 60):
    n = int(n0); us, ps = [], []
    while n > 1:
        us.append(np.log(n)); ps.append(n % 2)
        n = n // 2 if n % 2 == 0 else (3 * n + 1) // 2
    us = np.array(us); ps = np.array(ps)
    ax.scatter(us[:-1], np.diff(us), c=np.where(ps[:-1] == 1, PAL['coral'], PAL['teal']), s=2, lw=0, rasterized=True)
x = np.linspace(0.5, 16.5, 200)
ax.plot(x, np.log(1.5 + 0.5 * np.exp(-x)), color=PAL['coral'], lw=0.8)
ax.axhline(-np.log(2), color=PAL['teal'], lw=0.8)
ax.axhline(np.log(np.sqrt(3) / 2), color=PAL['ink'], lw=0.8, ls='--')
ax.text(1.0, 0.55, r'odd: $\log(3/2+1/(2x))$', ha='left', fontsize=7, color=PAL['coral'])
ax.text(1.0, -0.80, r'even: $-\log 2$', ha='left', fontsize=7, color=PAL['teal'])
ax.text(1.0, -0.10, r'mean drift $\delta^*$', ha='left', fontsize=7, color=PAL['ink'])
ax.set_xlabel(r'$u_n=\log x_n$'); ax.set_ylabel(r'$u_{n+1}-u_n$'); ax.set_ylim(-0.9, 0.75)
ax.set_title('(a)', loc='left', fontsize=8)
ax = axs[1]
n = 27; us = []; ps = []
while n > 1:
    us.append(np.log(n)); ps.append(n % 2); n = n // 2 if n % 2 == 0 else (3 * n + 1) // 2
us.append(0.0)
t = np.arange(len(us))
ax.plot(t, us, color=PAL['grey'], lw=0.7, zorder=1)
ax.scatter(t[:-1], us[:-1], c=np.where(np.array(ps) == 1, PAL['coral'], PAL['teal']), s=6, zorder=2, lw=0)
ax.plot(t, np.log(27) + t * np.log(np.sqrt(3) / 2), color=PAL['ink'], ls='--', lw=0.8)
ax.set_xlabel('step $n$ of $T$'); ax.set_ylabel(r'$u_n$'); ax.set_ylim(-0.5, 9.8)
ax.set_title('(b)', loc='left', fontsize=8)
fig.savefig(FIG + 'ch5b_switched.pdf'); plt.close(fig)

# 2. inverse tree counts and certified exponents
d = json.load(open(BOOK + '/scripts/ch5b_counts.json'))
Xs = np.array(d['Xs'], float); out = d['out']
fig, axs = plt.subplots(1, 2, figsize=(WIDTH, 2.2), gridspec_kw=dict(wspace=0.38))
ax = axs[0]
cols = plt.get_cmap(CMAP_SEQ)(np.linspace(0.1, 0.9, 6))
for i, root in enumerate(['4', '5', '7', '11', '17', '23']):
    rec = out[root]
    N = np.array(rec['N'], float); R = np.array(rec['R'], float)
    m = N > 0; ax.loglog(Xs[m], N[m], color=cols[i], lw=1.2, label=f'$n={root}$')
    m = R > 0; ax.loglog(Xs[m], R[m], color=cols[i], lw=1.0, ls='--')
xx = np.geomspace(1e2, 1e6, 10)
ax.loglog(xx, 0.6 * xx, color=PAL['grey'], lw=0.7, ls=':'); ax.text(1.5e3, 4e3, 'slope 1', fontsize=7, color=PAL['grey'])
ax.loglog(xx, 3.5 * xx ** 0.5, color=PAL['grey'], lw=0.7, ls=':'); ax.text(1.5e4, 120, 'slope 1/2', fontsize=7, color=PAL['grey'])
ax.set_xlabel(r'$X$'); ax.set_ylabel('count'); ax.legend(fontsize=6, ncol=1, loc='upper left')
ax.set_title('(a)', loc='left', fontsize=8)
ax = axs[1]
k = np.array([2, 3, 4, 5, 6, 7]); g = np.array([0.416, 0.435, 0.449, 0.460, 0.470, 0.477])
ax.plot(k, g, 'o-', color=PAL['violet'], label='doubling chain system')
ax.axhline(3 / 7, color=PAL['gold'], ls='--', lw=0.9); ax.text(4.6, 3/7 - 0.035, 'Krasikov, $3/7$', fontsize=7, color=PAL['gold'])
ax.axhline(0.84, color=PAL['coral'], ls='--', lw=0.9); ax.text(2.05, 0.848, 'Krasikov and Lagarias, $0.84$', fontsize=7, color=PAL['coral'])
ax.axhline(1.0, color=PAL['grey'], ls=':', lw=0.9); ax.text(2.05, 0.975, 'conjecture, $1$', fontsize=7, color=PAL['grey'])
ax.set_xlabel(r'modulus $3^k$, $k$'); ax.set_ylabel(r'certified exponent $\gamma_k$'); ax.set_ylim(0.38, 1.03)
ax.set_title('(b)', loc='left', fontsize=8)
fig.savefig(FIG + 'ch5b_tree.pdf'); plt.close(fig)

# 3. joint spectral radius
fig, axs = plt.subplots(1, 2, figsize=(WIDTH, 2.2), gridspec_kw=dict(wspace=0.35))
A = np.geomspace(1e-3, 0.995, 300)
for ax, b, ttl in ((axs[0], b0, r'(a) $\beta=0.585$'), (axs[1], 0.9, r'(b) $\beta=0.9$')):
    J = np.array([jsr_exact(a, b)[0] for a in A])
    rD = np.array([rho(mats(a, b)[0]) for a in A]); rN = np.array([rho(mats(a, b)[1]) for a in A])
    rP = np.array([np.sqrt(rho(mats(a, b)[1] @ mats(a, b)[0])) for a in A])
    ax.semilogx(A, J, color=PAL['plum'], lw=2.2, label='joint spectral radius')
    ax.semilogx(A, rD, color=PAL['green'], lw=1.5, ls=':', label='descent only (GD)')  # green and coral print as one gray: GD is dotted
    ax.semilogx(A, rN, color=PAL['sky'], lw=1.0, ls='--', label='injection only (NAG)')
    ax.semilogx(A, rP, color=PAL['coral'], lw=1.0, label='parity (SMM)')
    ax.set_xlabel(r'$a=\alpha\lambda$'); ax.set_ylim(0, 1.02); ax.set_title(ttl, loc='left', fontsize=8)
axs[0].set_ylabel('per step contraction'); axs[0].legend(fontsize=6, loc='lower left')
fig.savefig(FIG + 'ch5b_jsr.pdf'); plt.close(fig)

# 4. two step contraction on F_{mu,L}: bound, PEP and gradient descent
P = np.load(BOOK + '/scripts/ch5b_pep.npy')
fig, ax = plt.subplots(figsize=(WIDTH, 2.4))
K = np.geomspace(1.5, 200, 200)
for b, c, lab, ls, mk in ((0.3, PAL['sky'], r'$\beta=0.3$', (0, (4, 1.2, 1, 1.2)), 's'), (b0, PAL['coral'], r'$\beta=0.585$', '-', 'o'), (0.9, PAL['violet'], r'$\beta=0.9$', '-', 'o')):  # sky and coral print as one gray: beta=0.3 is dash-dot with squares
    q = (1 - 1 / K) * np.maximum(b, 1 - (1 + b) / K)
    ax.semilogx(K, -np.log(q) / 2 * K, color=c, lw=1.4, ls=ls, label=lab)
    rows = P[np.isclose(P[:, 0], b)]
    ax.semilogx(rows[:, 1], -np.log(np.sqrt(rows[:, 3])) / 2 * rows[:, 1], mk, color=c, ms=3.5 if mk == 's' else 4, mfc='white')
ax.semilogx(K, -np.log(1 - 1 / K) * K, color=PAL['green'], lw=1.2, ls='--', label='gradient descent')
ax.axhline((2 + b0) / 2, color=PAL['coral'], lw=0.6, ls=':')
ax.text(150, (2 + b0) / 2 + 0.03, r'$(2+\beta)/2$', fontsize=7, color=PAL['coral'], ha='right')
ax.set_xlabel(r'condition number $\kappa$'); ax.set_ylabel(r'rate per step $\times\ \kappa$')
ax.legend(fontsize=7, loc='lower right', ncol=2); ax.set_ylim(0, 1.9)
fig.savefig(FIG + 'ch5b_rate.pdf'); plt.close(fig)
print('done')
