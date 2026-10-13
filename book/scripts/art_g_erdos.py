"""Gallery: the Erdos 125 sumset.

A = numbers whose base 3 digits are all 0 or 1, B = numbers whose base 4 digits
are all 0 or 1, and r(n) = #{(a, b) in A x B : a + b = n}.
Top: r(n) for 0 <= n < 3^12 = 531441 in a base 3 Morton layout on a 729 x 729
square: if n has base 3 digits d_0 d_1 ... d_11 (d_0 least significant), it
is drawn in column sum d_{2i} 3^i and row sum d_{2i+1} 3^i (row 0 at the
bottom).  Colour is log r(n); background where n is not in A + B.  In this
layout A itself is the product Cantor dust {0,1}-digit set in both axes.
Bottom: the density |(A+B) cap [0, N]| / (N+1) for 10^3 <= N <= 4^11, with the
powers of 3 marked and the Hasler and Melfi bound 1015/1458 (dashed).
"""
from art_common import *
from art_coral import EMBER
from matplotlib.colors import ListedColormap, BoundaryNorm

E = 11
NN = 4 ** E
def digitset(base, limit):
    s = np.array([0], dtype=np.int64); p = 1
    while p <= limit:
        s = np.concatenate([s, s + p]); p *= base
    return np.sort(s[s <= limit])
A = digitset(3, NN); B = digitset(4, NN)
rep = np.zeros(NN + 1, dtype=np.int32)
for b in B:
    a = A[A <= NN - b]; rep[a + b] += 1
mark = rep > 0
cnt = np.cumsum(mark, dtype=np.int64)
dens = lambda n: cnt[n] / (n + 1)
print('|A|,|B| up to 4^11:', len(A), len(B))
for k in (8, 9, 10, 11): print('density at 4^%d: %.4f' % (k, dens(4 ** k)))
for m in range(7, 13): print('density at 3^%d: %.4f   at 3^%d - 1: %.4f' % (m, dens(3 ** m), m, dens(3 ** m - 1)))
lo = 10 ** 3
nn = np.arange(lo, NN + 1); d = cnt[lo:] / (nn + 1.0); i = int(np.argmin(d))
print('min density over [1e3, 4^11]: %.4f at N = %d' % (d[i], lo + i))
M = 3 ** 12
print('max r(n) for n < 3^12:', rep[:M].max(), 'at', np.argmax(rep[:M]), ' fraction of n < 3^12 in A+B: %.4f' % mark[:M].mean())
print('r distribution:', np.bincount(rep[:M])[:8] / M)

Wd = 729
n_ = np.arange(M); cx = np.zeros(M, np.int64); cy = np.zeros(M, np.int64); t = n_.copy()
for i in range(6):
    cx += (t % 3) * 3 ** i; t //= 3
    cy += (t % 3) * 3 ** i; t //= 3
img = np.zeros((729, 729)); img[cy, cx] = rep[:M]
CM = mpl.colors.LinearSegmentedColormap.from_list('e2', EMBER(np.linspace(0.08, 1, 256)))
CM.set_bad(NIGHT)
rmax = rep[:M].max()
show = np.where(img > 0, img, np.nan)
lnorm = mpl.colors.LogNorm(1, rmax)
fig = plt.figure(figsize=(4.5, 6.6), facecolor=NIGHT)
side = 0.86
ax0 = fig.add_axes([0.04, 0.395, side, side * 4.5 / 6.6])
ax0.imshow(show, cmap=CM, norm=lnorm, interpolation='nearest', aspect='auto', origin='lower')
ax0.axis('off')
cax = fig.add_axes([0.915, 0.395, 0.022, side * 4.5 / 6.6])
cb = fig.colorbar(mpl.cm.ScalarMappable(norm=lnorm, cmap=CM), cax=cax)
cb.outline.set_visible(False); cb.ax.tick_params(colors='#C9CCD8', labelsize=7, which='both')
cb.set_ticks([1, 2, 5, 10, 20, 50]); cb.set_ticklabels(['1', '2', '5', '10', '20', '50'])
cb.ax.minorticks_off()
cb.set_label(r'$r(n)$', color='#E6E6EE')
ax = fig.add_axes([0.13, 0.065, 0.83, 0.22]); ax.set_facecolor(NIGHT); ax.grid(False)
for sp in ('left', 'bottom'):
    ax.spines[sp].set_color('#8D99AE')
ax.tick_params(colors='#C9CCD8')
g = np.unique(np.round(np.logspace(3, np.log10(NN), 5000)).astype(np.int64))
g = np.unique(np.concatenate([g, [3 ** m - 1 for m in range(7, 15) if 3 ** m - 1 <= NN]]))
for m in range(7, 15):
    if 3 ** m <= NN:
        ax.axvline(3 ** m, color=PAL['gold'], lw=0.5, alpha=0.5)
        ax.text(3 ** m, 1.005, f'$3^{{{m}}}$', color=PAL['gold'], fontsize=6.5, ha='center', va='bottom')
from matplotlib.collections import LineCollection
pts = np.c_[np.log10(g), dens(g)]
segs = np.stack([pts[:-1], pts[1:]], 1)
segs[:, :, 0] = 10 ** segs[:, :, 0]
lc = LineCollection(segs, cmap=EMBER, norm=mpl.colors.Normalize(0.6, 1.0), lw=1.1)
lc.set_array(dens(g[:-1])); ax.add_collection(lc)
ax.axhline(1015 / 1458, color='#7FE0FF', lw=0.7, ls=(0, (4, 3)))
ax.text(1.15e3, 1015 / 1458 - 0.012, 'Hasler and Melfi bound $1015/1458$ on the lower density', color='#7FE0FF', fontsize=6.5, va='top')
ax.set_xscale('log'); ax.set_xlim(1e3, NN); ax.set_ylim(0.6, 1.0)
ax.set_xlabel(r'$N$', color='#E6E6EE')
ax.set_ylabel(r'density of $A+B$ in $[0,N]$', color='#E6E6EE', fontsize=8)
save(fig, 'art_g_erdos', facecolor=NIGHT, bbox_inches=None)
