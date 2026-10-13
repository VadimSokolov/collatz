import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
import sys; sys.path.insert(0, BOOK + '/scripts')
from bookstyle import *
import numpy as np, math, json
from matplotlib.patches import FancyArrowPatch
FIG = BOOK + '/figures/'
L32, L2, L3 = math.log(1.5), math.log(2), math.log(3)
def T(n): return n // 2 if n % 2 == 0 else (3 * n + 1) // 2
# (a) residue graph mod 32 on a circle; edges r -> lifts; three of its 91 positive simple cycles highlighted (those of the negative integer cycles)
k = 5; M = 2 ** k
fig, axs = plt.subplots(1, 2, figsize=(WIDTH, 2.4), gridspec_kw=dict(width_ratios=[1.25, 1], wspace=0.32))
ax = axs[0]
ang = {r: 2 * np.pi * r / M for r in range(M)}
pos = {r: (np.cos(ang[r]), np.sin(ang[r])) for r in range(M)}
for r in range(M):
    t = T(r) % (M // 2)
    for lift in (t, t + M // 2):
        ax.annotate('', xy=pos[lift], xytext=pos[r], arrowprops=dict(arrowstyle='-', color='#D8D8D8', lw=0.5,
                    connectionstyle='arc3,rad=0.25'))
cyc_cols = [PAL['coral'], PAL['gold'], PAL['violet']]
for cyc, col, lab in (([-1], PAL['coral'], r'$-1$: mean $\log(3/2)$'), ([-5, -7, -10], PAL['gold'], r'$-5$ cycle: mean $0.039$'),
                      ([-17, -25, -37, -55, -82, -41, -61, -91, -136, -68, -34], PAL['violet'], r'$-17$ cycle: mean $0.006$')):
    rs = [c % M for c in cyc]
    for i in range(len(rs)):
        a, b = rs[i], rs[(i + 1) % len(rs)]
        if a == b:
            ax.plot(pos[a][0] * 1.12, pos[a][1] * 1.12, 'o', color=col, ms=9, mfc='none', mew=1.6)
        else:
            ax.annotate('', xy=pos[b], xytext=pos[a], arrowprops=dict(arrowstyle='-|>', color=col, lw=1.4, mutation_scale=8,
                        connectionstyle='arc3,rad=0.25'))
    ax.plot([], [], color=col, lw=1.6, label=lab)
for r in range(M):
    ax.plot(*pos[r], 'o', color=PAL['teal'] if r % 2 == 0 else PAL['coral'], ms=3.2, zorder=5)
    ax.text(pos[r][0] * 1.17, pos[r][1] * 1.17, str(r), fontsize=5, ha='center', va='center', color=PAL['ink'])
ax.set_aspect('equal'); ax.axis('off'); ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.3, 1.3)
ax.legend(fontsize=6, loc='upper center', bbox_to_anchor=(0.5, 0.02), ncol=1, handlelength=1.4)   # below the circle, clear of the labels
ax.set_title('(a) residue graph modulo $32$', loc='left', fontsize=8)
# (b) block descent frequencies and expected block change
ax = axs[1]
d = json.load(open(BOOK + '/scripts/data/ch5c_blocks.json'))
js = np.array(sorted(int(j) for j in d)); fr = np.array([d[str(j)][0] for j in js]); mod = np.array([d[str(j)][1] for j in js])
ax.bar(js - 0.2, fr, 0.4, color=PAL['teal'], label='observed, $q<2^{16}$')
ax.bar(js + 0.2, mod, 0.4, color=PAL['gold'], label=r'$2^{-\lceil j\log_2(3/2)\rceil+1}$')
ax.set_yscale('log'); ax.set_xlabel('run length $j$ (ascents $j-1$)'); ax.set_ylabel('fraction of blocks that descend')
ax2 = ax.twinx(); ax2.plot(js, (js - 1) * L32 + L3 - 3 * L2, 's-', color=PAL['coral'], ms=3, lw=0.9, label='mean block change')
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, fontsize=6, loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=1)   # below the axes, clear of the markers
ax2.axhline(0, color=PAL['coral'], lw=0.5, ls=':'); ax2.set_ylabel('mean block change', color=PAL['coral'])
ax2.tick_params(axis='y', colors=PAL['coral']); ax2.grid(False); ax2.spines['right'].set_visible(True)
ax.set_title('(b) blocks', loc='left', fontsize=8)
fig.savefig(FIG + 'ch5c_residue_blocks.pdf'); plt.close(fig)
# (c) failure fraction of the run length potential as a function of c (model), with scan points
def pmodel(c):
    p = 0.0
    for kk in range(2, 80):
        need = (kk * L2 - L3) / c; t = math.ceil(need)     # a failure needs v2(y+1) - 1 >= need
        p += 2.0 ** (-kk) * 2.0 ** (-t)
    return p
cs = np.linspace(0.41, 3.0, 400); pm = np.array([pmodel(c) for c in cs])
fig, ax = plt.subplots(figsize=(WIDTH, 2.1))
ax.plot(cs, pm, color=PAL['violet'], lw=1.4, label='model')
ax.plot([0.5, 1.0, 2.0], [0.1613, 0.2087, 0.2333], 'o', color=PAL['coral'], ms=5, label=r'scan of odd $x<2^{20}$')
ax.axvline(L32, color=PAL['grey'], lw=0.7, ls='--'); ax.text(L32 + 0.03, 0.05, r'$c=\log(3/2)$', fontsize=7, color=PAL['grey'])
ax.set_xlabel('weight $c$'); ax.set_ylabel('fraction of odd $x$\n' + r'where $V_c$ fails'); ax.legend(fontsize=7, loc='lower right')
ax.set_ylim(0, 0.3)
fig.savefig(FIG + 'ch5c_runpotential.pdf'); plt.close(fig)
print('inf over c of model failure:', pm.min(), 'at c =', cs[pm.argmin()])
