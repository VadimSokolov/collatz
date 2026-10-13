import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
# Figure for the gain weight theorem (Chapter 7): (a) exact matched speed floor ratio against a, converging to
# F = p sum c^2 / M^2; (b) slow mode speed M/p against F for momentum on every proper subset of residues, p <= 6.
import sys; sys.path.insert(0, BOOK + '/scripts')
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    import ch7_gainweight as gw
from bookstyle import *
import numpy as np, json

b0 = gw.b0
sel = [('GD', [gw.D], PAL['green']),
       (r'every step NAG, $\beta=0.9$', [gw.NAG(0.9)], PAL['sky']),
       (r'SMM, $\beta=0.585$', [gw.D, gw.NAG(b0)], PAL['gold']),
       (r'SMM, $\beta=0.9$', [gw.D, gw.NAG(0.9)], PAL['coral']),
       (r'$p=3$, HB on $\{1,2\}$, $\beta=0.9$', [gw.D, gw.HB(0.9), gw.HB(0.9)], PAL['violet']),
       (r'cyclic steps $(1,3)$', [gw.GDh(1.0), gw.GDh(3.0)], PAL['plum']),
       (r'generic, $p=3$, $K=3$', gw.gen, PAL['teal'])]
A = np.logspace(-4, -0.7, 60)
fig, axs = plt.subplots(2, 1, figsize=(WIDTH, 5.4), gridspec_kw=dict(height_ratios=[1, 1], hspace=0.75))
ax = axs[0]
# Grayscale: coral, green and sky print alike, as do teal, violet and plum, so every curve gets its own line style and marker
# (the limits F are drawn dotted and no curve is dotted). Colours are unchanged.
LS = {'solid': '-', 'dash': (0, (3.2, 1.4)), 'dashdot': (0, (3.2, 1.1, 0.7, 1.1))}
sty = [('solid', None), ('dash', 'v'), ('solid', 'o'), ('dashdot', 's'), ('dash', 'D'), ('dashdot', '^'), ('solid', 'X')]
MS = {None: 0, 'o': 4.4, 's': 4.0, 'D': 4.0, '^': 4.8, 'v': 4.8, 'X': 4.8}
hl = []
for k, ((name, sch, c), (ls, mkr)) in enumerate(zip(sel, sty)):
    p = len(sch); cs = gw.gains_by_impulse(sch, 3000); M = cs.sum(); F = p * (cs ** 2).sum() / M ** 2
    vals = []
    for a in A:
        v, r = gw.exact_variance(sch, a); hg = M / p
        vals.append(v / (hg / (a * (2 - a * hg))) if np.isfinite(v) else np.nan)
    ax.axhline(F, color=c, lw=1.0, ls=(0, (1, 2)), alpha=0.9)
    kw = dict(color=c, lw=1.9, ls=LS[ls], marker=mkr, ms=MS[mkr], mec='white', mew=0.5)
    ax.plot(A, vals, markevery=(2 + k, 9), **kw)
    hl.append(plt.Line2D([], [], label=name, **kw))
ax.set_xscale('log'); ax.set_xlim(1e-4, 0.2); ax.set_ylim(0.85, 1.38)
ax.set_xlabel(r'normalized step $a=\alpha\lambda$'); ax.set_ylabel('floor relative to GD\nat matched speed')
ax.legend(handles=hl, loc='upper center', bbox_to_anchor=(0.5, -0.24), fontsize=7, ncol=3, columnspacing=0.7, handlelength=2.8, handletextpad=0.5)
ax.text(-0.15, 1.04, '(a)', transform=ax.transAxes, va='bottom')
ax = axs[1]
d = json.load(open(BOOK + '/scripts/ch7_gainweight.json'))
pcol = {2: PAL['gold'], 3: PAL['coral'], 4: PAL['plum'], 5: PAL['violet'], 6: PAL['sky']}
mk = {0.5: 'o', 0.585: 's', 0.9: '^'}
# Grayscale: p = 3 and 6 (coral, sky) and p = 4 and 5 (plum, violet) print alike, so the marker size grows with p and p = 4, 6 are drawn open.
psz = {2: 10, 3: 14, 4: 20, 5: 26, 6: 34}
pfill = {2: True, 3: True, 4: False, 5: True, 6: False}
for (b, p, mask, spd, F) in sorted(d['subset_points'], key=lambda q: -q[1]):
    if sum(mask) in (0, p): continue
    if pfill[p]: ax.scatter(spd, F, s=psz[p], marker=mk[round(b, 3)], color=pcol[p], alpha=0.8, lw=0, zorder=3)
    else: ax.scatter(spd, F, s=psz[p], marker=mk[round(b, 3)], facecolors='none', edgecolors=pcol[p], alpha=0.9, linewidths=0.9, zorder=3)
for b in [0.5, b0, 0.9]:
    ax.scatter([1 / (1 - b)], [1.0], s=34, marker=mk[round(b, 3)], facecolor='white', edgecolor=PAL['green'], lw=1.3, zorder=4)
ax.scatter([1.0], [1.0], s=34, marker='D', color=PAL['green'], zorder=4)
ax.axhline(1.0, color=PAL['green'], lw=1.0)
ax.text(2.6, 0.975, r'GD and every step momentum, $F=1$', fontsize=7, color=PAL['green'], va='top')
ax.set_xscale('log'); ax.set_xlabel(r'slow mode speed multiplier $M/p$')
ax.set_ylabel(r'floor factor $F$'); ax.set_ylim(0.93, 1.42); ax.set_xlim(0.9, 12)
ax.set_xticks([1, 2, 5, 10]); ax.set_xticklabels(['1', '2', '5', '10']); ax.minorticks_off()
hp = [plt.Line2D([], [], marker='o', ls='', color=pcol[p], ms=psz[p] ** 0.5, mfc=pcol[p] if pfill[p] else 'none', mew=0 if pfill[p] else 0.9, label=rf'$p={p}$') for p in range(2, 7)]
hb = [plt.Line2D([], [], marker=m, ls='', color=PAL['grey'], ms=4.5, label=rf'$\beta={bb}$') for bb, m in [('0.5', 'o'), ('0.585', 's'), ('0.9', '^')]]
l1 = ax.legend(handles=hp, loc='upper right', fontsize=7, ncol=5, columnspacing=0.6, handletextpad=0.1)
ax.add_artist(l1); ax.legend(handles=hb, loc='upper right', bbox_to_anchor=(1.0, 0.86), fontsize=7, ncol=3, columnspacing=0.6, handletextpad=0.1)
ax.text(-0.15, 1.04, '(b)', transform=ax.transAxes, va='bottom')
fig.savefig(BOOK + '/figures/ch7_gainweight.pdf')
