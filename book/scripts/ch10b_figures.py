import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
import sys; sys.path.insert(0, BOOK + '/scripts')
from bookstyle import *
import numpy as np, json, math
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
FIG = BOOK + '/figures/'
d = json.load(open(BOOK + '/scripts/data/ch10b_sieve.json'))
k = np.array([r[0] for r in d]); cnt = np.array([float(r[1]) for r in d]); prop = np.array([r[2] for r in d])
H = lambda p: -p * np.log2(p) - (1 - p) * np.log2(1 - p)
h = H(1 / math.log2(3))
fig, axs = plt.subplots(1, 2, figsize=(WIDTH, 2.3), gridspec_kw=dict(wspace=0.42))
ax = axs[0]
ax.semilogy(k, cnt, color=PAL['coral'], lw=1.6, label='uncovered classes')
ax.semilogy(k, 2.0 ** k, color=PAL['grey'], lw=0.8, ls=':', label='all classes $2^k$')
ax.semilogy(k, 2.0 ** (h * k) * cnt[-1] / 2.0 ** (h * k[-1]), color=PAL['violet'], lw=0.9, ls='--', label=r'$2^{0.95k}$')
ax.set_xlabel('level $k$'); ax.set_ylabel('count'); ax.legend(fontsize=6.5, loc='upper left')
ax.set_ylim(0.5, 1e125); ax.set_yticks([1, 1e30, 1e60, 1e90, 1e120])   # counts are at least 1
ax.set_title('(a)', loc='left', fontsize=8)
# (b) icicle of the 2 adic exceptional set approximations, levels 1..16
ax = axs[1]
K = 16; M = 2 ** K
def iota(r, kk):
    s = 0.0
    for i in range(kk): s += ((r >> i) & 1) * 2.0 ** (-i - 1)
    return s
alive_level = np.full(M, K + 1)      # first level at which the class stops (coefficient stopping time)
for r in range(M):
    x = r; m = 0
    for j in range(1, K + 1):
        if x & 1: m += 1; x = (3 * x + 1) >> 1
        else: x >>= 1
        if 3 ** m < 2 ** j: alive_level[r] = j; break
cmap = plt.get_cmap('viridis')
patches = []; cols = []
for kk in range(1, K + 1):
    mod = 2 ** kk
    seen = set()
    for r in range(M):
        if alive_level[r] > kk:
            rr = r % mod
            if rr in seen: continue
            seen.add(rr)
            patches.append(Rectangle((iota(rr, kk), kk - 0.46), 2.0 ** (-kk), 0.92))   # level kk drawn at height kk
            cols.append(cmap(0.05 + 0.9 * (kk - 1) / (K - 1)))
pc = PatchCollection(patches, facecolors=cols, edgecolors='none'); ax.add_collection(pc)
ax.set_xlim(0, 1); ax.set_ylim(K + 0.5, 0.5); ax.set_yticks([1, 4, 8, 12, 16])
ax.set_xlabel(r'$2$ adic coordinate $\iota(r)$'); ax.set_ylabel('level $k$')
ax.axvline(1.0 - 2 ** -K, color=PAL['sky'], lw=0.6)
ax.text(1.0, 1.02, r'$-1$', transform=ax.transAxes, color=PAL['sky'], fontsize=7, ha='center', va='bottom')
ax.grid(False); ax.set_title('(b)', loc='left', fontsize=8)
fig.savefig(FIG + 'ch10b_sieve.pdf'); plt.close(fig)
print('h =', h)
