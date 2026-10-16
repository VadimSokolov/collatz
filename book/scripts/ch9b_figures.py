"""Figures and renewal constants for the chapter on positive density from the backward tree (adapted from N. Polson's
renewal_and_fourier_fig.py, collatz_steps.py, tree.py and fc.py). Reads data/ch9b_steps.json and data/ch9b_offset_law.json.
ch9b_steps.pdf   : histograms of sigma0(n)/log n and A(n)/log n for 10^6 <= n < 3*10^6, with the model constants and Lau's.
ch9b_fourier.pdf : (a) the largest Fourier coefficient of the offset law at frequencies prime to 3, against d^-2 and a
                   geometric curve; (b) the distances of rho_d from the lifts of rho_{d-1} and of rho_{ceil(d/2)}.
ch9b_tree.pdf    : the odd integers below 400 within five Syracuse steps of 1, by the length of their word and log_2 n.
ch9b_crossing.pdf: fourteen random words drawn (letters Geometric(1/2), seed 3), those crossing within 80 letters shown, and their first crossings of the line
                   of slope log_2 3 through (delta, 2 delta), delta = 40, after length 2 delta / 5.
Prints the profile |mu_d^(2^k)| for d >= 14 and the renewal constants: the standard deviations of A/log n and sigma0/log n times sqrt(log n), the margin of 10.46
over 3/log(4/3), and the log n and the number of decimal digits at which the margin covers 90 and 99 percent."""
import sys, os, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bookstyle import *
import numpy as np
from statistics import NormalDist
HERE = os.path.dirname(os.path.abspath(__file__)); FIG = os.path.join(HERE, '..', 'figures'); DATA = os.path.join(HERE, 'data')

# ---------------------------------------------------------------- renewal constants
mu = math.log(4 / 3); s2 = 2 * math.log(2)**2; k = 1 + math.log2(3)
sdA = math.sqrt(s2 / mu**3); sdT = k * sdA; gap = 10.46 - 3 / mu
print(f'sd(A/L)*sqrt(L) = {sdA:.4f}, sd(sigma0/L)*sqrt(L) = {sdT:.4f}, 3/log(4/3) = {3 / mu:.5f}, 1/log(4/3) = {1 / mu:.5f}, margin = {gap:.5f}')
print(f'(1+log2 3)*3.4881 + 1/log 2 = {k * 3.4881 + 1 / math.log(2):.5f}')
for p in (0.9, 0.99):
    L = (NormalDist().inv_cdf(p) * sdT / gap)**2
    print(f'coverage {p}: log n = {L:.4g}, digits = {L / math.log(10):.4g}')
st = json.load(open(os.path.join(DATA, 'ch9b_steps.json')))
for key, r in st['ranges'].items():
    print(key, {q: round(v, 4) if isinstance(v, float) else v for q, v in r.items()}, 'predicted sd', round(sdT / math.sqrt(r['logmid']), 3))

# ---------------------------------------------------------------- step histograms
fig, axs = plt.subplots(1, 2, figsize=(WIDTH, 1.9), gridspec_kw=dict(wspace=0.32))
for ax, key, col, model, lau, lab, title in [(axs[0], 'T', PAL['teal'], 3 / mu, 10.46, r'$\sigma_0(n)/\log n$', '(a) steps of $T_0$'),
                                              (axs[1], 'A', PAL['green'], 1 / mu, 3.4881, r'$A(n)/\log n$', '(b) steps of $S$')]:
    h = st['hist'][key]; e = np.array(h['edges'])
    ax.stairs(h['dens'], e, fill=True, color=col, alpha=0.75, lw=0)
    ax.axvline(model, color=PAL['coral'], lw=0.9, ls='--', label=f'model {model:.3f}')
    ax.axvline(lau, color=PAL['ink'], lw=0.9, label=f'Lau {lau}')
    ax.set_xlabel(lab); ax.set_title(title, loc='left', fontsize=8); ax.legend(fontsize=5.5, loc='upper right')
axs[0].set_ylabel('density')
fig.savefig(os.path.join(FIG, 'ch9b_steps.pdf')); plt.close(fig)

# ---------------------------------------------------------------- offset law
r = json.load(open(os.path.join(DATA, 'ch9b_offset_law.json')))
d = np.array([x['d'] for x in r]); f = np.array([x['fmax'] for x in r])
print('d^2 fmax', np.round(f * d**2, 3)); print('ratios', np.round(f[1:] / f[:-1], 3))
print('worst xi a power of 2 at every level:', all(x['pow2'] for x in r), [round(x['log2xi'], 2) for x in r])
print('worst k / d:', [round(x['log2xi'] / x['d'], 3) for x in r if x['d'] >= 8], 'log_3(3^d / xi):', [round(x['d'] - x['log2xi'] * math.log(2) / math.log(3), 2) for x in r if x['d'] >= 8])
for x in r:
    if 'profile' in x:
        p = {int(j): v for j, v in x['profile'].items()}; top = max(p.values())
        print(f"d={x['d']} profile k=16..27:", ' '.join(f'{p[j]:.4f}' for j in range(16, 28) if j in p), ' k with at least half the peak:', [j for j in sorted(p) if p[j] >= top / 2])
fig, axs = plt.subplots(1, 2, figsize=(WIDTH, 2.0), gridspec_kw=dict(wspace=0.36))
ax = axs[0]
ax.loglog(d, f, 'o-', color=PAL['teal'], ms=2.5, lw=0.9, label=r'$\max_{3\nmid\xi}|\widehat\pi_d(\xi)|$')
ax.loglog(d, 3.7 * d**-2.0, '--', color=PAL['coral'], lw=0.9, label=r'$3.7\,d^{-2}$')
ax.loglog(d, 0.58 * 0.8**(d - 1), ':', color=PAL['grey'], lw=1.1, label=r'$0.58\cdot0.8^{d-1}$')
ax.set_xlabel('level $d$'); ax.set_title('(a)', loc='left', fontsize=8); ax.legend(fontsize=5.5, loc='lower left')
ax = axs[1]
ax.loglog(d, [x['L1'] for x in r], 'o-', color=PAL['green'], ms=2.5, lw=0.9, label=r'from $\rho_{d-1}$')
ax.loglog(d[1:], [x['L1half'] for x in r][1:], 's--', color=PAL['plum'], ms=2.5, lw=0.9, label=r'from $\rho_{\lceil d/2\rceil}$')
ax.set_xlabel('level $d$'); ax.set_ylabel('mean distance'); ax.set_title('(b)', loc='left', fontsize=8); ax.legend(fontsize=5.5, loc='lower left')
fig.savefig(os.path.join(FIG, 'ch9b_fourier.pdf')); plt.close(fig)

# ---------------------------------------------------------------- the backward tree from 1
def syr(y):
    y = 3 * y + 1
    return y >> ((y & -y).bit_length() - 1)
depth = {1: 0}
def dep(n):
    if n not in depth: depth[n] = 1 + dep(syr(n))
    return depth[n]
small = [n for n in range(1, 400, 2) if dep(n) <= 5]
print('tree figure:', len(small), 'nodes, depths', sorted({depth[n] for n in small}))
fig, ax = plt.subplots(figsize=(WIDTH, 2.7))
for n in small[1:]:
    ax.plot([depth[n] - 1, depth[n]], [math.log2(syr(n)), math.log2(n)], color=PAL['grey'], lw=0.5, zorder=1)
ax.scatter([depth[n] for n in small], [math.log2(n) for n in small], s=9, color=PAL['teal'], zorder=2)
for n in small:
    if n <= 100: ax.annotate(str(n), (depth[n], math.log2(n)), fontsize=5.5, xytext=(3, 0), textcoords='offset points', va='center')
ax.set_xlabel('Syracuse steps to 1 (length of the word)'); ax.set_ylabel(r'$\log_2 n$')
fig.savefig(os.path.join(FIG, 'ch9b_tree.pdf')); plt.close(fig)

# ---------------------------------------------------------------- first crossings
rng = np.random.default_rng(3); delta = 40; s3 = math.log2(3); n = np.arange(1, 81)
fig, ax = plt.subplots(figsize=(WIDTH, 2.7))
ax.plot(np.arange(0, 80), 2 * delta + (np.arange(0, 80) - delta) * s3, color=PAL['coral'], lw=1.0, label=r'slope $\log_2 3$ through $(\delta,2\delta)$')
lens = []
for _ in range(14):
    deg = np.cumsum(rng.geometric(0.5, size=80))
    hits = np.where((deg >= 2 * delta + (n - delta) * s3) & (n > delta - 3 * delta // 5))[0]
    if len(hits) == 0: continue
    c = int(hits[0]) + 1; lens.append(c)
    ax.plot(n[:c], deg[:c], color=PAL['teal'], lw=0.7, alpha=0.85)
    ax.scatter([c], [deg[c - 1]], color=PAL['ink'], s=8, zorder=3)
ax.axvline(delta, color=PAL['grey'], lw=0.5, ls=':'); ax.axhline(2 * delta, color=PAL['grey'], lw=0.5, ls=':')
ax.set_xlim(0, 80); ax.set_ylim(0, 150); ax.set_xlabel('length of the word'); ax.set_ylabel('degree (sum of the letters)')
ax.legend(fontsize=6, loc='upper left')
fig.savefig(os.path.join(FIG, 'ch9b_crossing.pdf')); plt.close(fig)
print('first crossings at delta = 40:', sorted(lens))
