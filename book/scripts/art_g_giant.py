"""Gallery: the Syracuse orbit of x0 = 3^3155, a 5001 bit number close to 2^5000.

y_j = S^(j)(x0), k_j = v2(3 y_{j-1} + 1), computed in exact integer arithmetic.
Top: the exponents k_1, ..., k_N as a mosaic read left to right, top to bottom
(132 per row).  Middle: log2 y_j, with the model drift line log2 x0 - (2 - log2 3) j
(dashed), which reaches 0 at j = 12050.  Bottom: the running mean
kbar_N = (k_1 + ... + k_N)/N with the levels log2 3 and 2.
"""
from art_common import *
from art_coral import EMBER
import math
from matplotlib.colors import ListedColormap, BoundaryNorm

x = 3 ** 3155
ks, lg = [], []
y = x
while y != 1:
    z = 3 * y + 1
    k = (z & -z).bit_length() - 1
    y = z >> k
    ks.append(k)
    lg.append(math.log2(y) if y < 2 ** 1000 else (y.bit_length() - 60) + math.log2(y >> (y.bit_length() - 60)))
import math as _m
ks = np.array(ks); lg = np.array(lg)
N = len(ks)
lx0 = (x.bit_length() - 60) + _m.log2(x >> (x.bit_length() - 60))
kbar = np.cumsum(ks) / np.arange(1, N + 1)
print('log2 x0 = %.4f' % lx0, 'strides', N, 'T0 steps', N + ks.sum(), 'T steps', ks.sum())
print('kbar_N final %.5f' % kbar[-1], ' predicted (log2 x0)/N + log2 3 ~ %.5f' % (lx0 / N + np.log2(3)))
print('frequency of k = 1..6:', np.round(np.bincount(ks, minlength=8)[1:7] / N, 4), 'max k', ks.max())
print('model line hits zero at j = %.1f' % (lx0 / (2 - np.log2(3))))
print('mean descent per stride %.4f (2 - log2 3 = %.4f)' % (lx0 / N, 2 - np.log2(3)))

W = 132
rows = int(np.ceil(N / W))
grid = np.full(rows * W, np.nan); grid[:N] = np.minimum(ks, 6); grid = grid.reshape(rows, W)
cols = [EMBER(c) for c in (0.12, 0.36, 0.55, 0.72, 0.87, 1.0)]
cmap = ListedColormap(cols); cmap.set_bad(NIGHT)
bn = BoundaryNorm(np.arange(0.5, 7.5), 6)

fig = plt.figure(figsize=(4.5, 6.6), facecolor=NIGHT)
def style(ax):
    ax.set_facecolor(NIGHT); ax.grid(False)
    for sp in ('left', 'bottom'):
        ax.spines[sp].set_color('#8D99AE')
    ax.tick_params(colors='#C9CCD8')
ax0 = fig.add_axes([0.03, 0.43, 0.94, 0.94 * W and 0.94 * 4.5 * rows / W / 6.6])
ax0.imshow(grid, cmap=cmap, norm=bn, interpolation='nearest', aspect='auto')
ax0.axis('off')
top = 0.43 + 0.94 * 4.5 * rows / W / 6.6
# legend chips
for i in range(6):
    xx = 0.22 + i * 0.11
    fig.patches.append(mpl.patches.Rectangle((xx, top + 0.012), 0.022, 0.014, transform=fig.transFigure,
                                             facecolor=cols[i], edgecolor='none'))
    fig.text(xx + 0.028, top + 0.019, f'$k={i + 1}$' if i < 5 else r'$k\geq 6$', color='#E6E6EE', fontsize=7, va='center')
ax1 = fig.add_axes([0.13, 0.235, 0.82, 0.15]); style(ax1)
j = np.arange(1, N + 1)
ax1.plot(np.r_[0, j], np.r_[lx0, lg], color=PAL['coral'], lw=0.8)
ax1.set_xlim(0, N); ax1.set_ylim(0, 5100)
ax1.set_ylabel(r'$\log_2 y_j$', color='#E6E6EE'); ax1.set_xticklabels([])
jm = lx0 / (2 - np.log2(3))
ax1.plot([0, jm], [lx0, 0], color='#FFF3C4', lw=0.5, ls=(0, (3, 3)), alpha=0.8)
ax2 = fig.add_axes([0.13, 0.065, 0.82, 0.15]); style(ax2)
ax2.plot(j, kbar, color='#7FE0FF', lw=0.9)
ax2.axhline(np.log2(3), color=PAL['gold'], lw=0.7)
ax2.axhline(2, color='#C9CCD8', lw=0.5, ls=(0, (3, 3)))
ax2.set_xlim(0, N); ax2.set_ylim(1.45, 2.25)
ax2.set_xlabel(r'stride $j$', color='#E6E6EE')
ax2.set_ylabel(r'$\bar k_j$', color='#E6E6EE')
ax2.text(N * 0.985, np.log2(3) - 0.03, r'$\log_2 3$', color=PAL['gold'], fontsize=7, ha='right', va='top')
ax2.text(N * 0.985, 2.03, r'$2$', color='#C9CCD8', fontsize=7, ha='right', va='bottom')
save(fig, 'art_g_giant', facecolor=NIGHT, bbox_inches=None)
