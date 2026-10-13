"""Gallery: the orbits of 1 <= n <= 2000 as a heat map.

Row n is the T_0 orbit of n (x -> 3x+1 odd, x/2 even), written right aligned
so that every orbit ends at 1 in the last column; the colour of a cell is
log2 of the iterate.  Cells to the left of the starting value are empty.  Rows are sorted by
orbit length (ties broken by n), so the left edge of the coloured region is
the sorted total stopping time and shared orbit tails line up as columns.
"""
from art_common import *
from art_coral import EMBER

N = 2000
orbs = [collatz_orbit(n) for n in range(1, N + 1)]
L = max(len(o) for o in orbs)
img = np.full((N, L), np.nan)
for i, o in enumerate(orbs):
    img[i, L - len(o):] = np.log2(np.array(o, float))
lens = np.array([len(o) - 1 for o in orbs])
rank = np.lexsort((np.arange(N), lens))
img = img[rank]
k = np.argmax(lens); print('longest orbit: n =', k + 1, 'steps', lens[k], 'columns', L)
mxv = max(max(o) for o in orbs); im = np.argmax([max(o) for o in orbs])
print('largest value', mxv, 'log2', np.log2(mxv), 'from n =', im + 1)
# fraction of rows that pass through 9232
print('orbits through 9232:', sum(9232 in o for o in orbs), 'steps 9232 to 1:', len(collatz_orbit(9232)) - 1)
print('lengths at quartile rows', np.sort(lens)[[0, 500, 1000, 1500, 1999]])

fig = plt.figure(figsize=(4.5, 5.6), facecolor=NIGHT)
ax = fig.add_axes([0.08, 0.1, 0.89, 0.8]); ax.set_facecolor(NIGHT)
norm = mpl.colors.Normalize(0, np.ceil(np.log2(mxv)))
ax.imshow(img, aspect='auto', cmap=EMBER, norm=norm, interpolation='nearest',
          extent=(-L + 0.5, 0.5, N + 0.5, 0.5))
ax.grid(False)
for sp in ('left', 'bottom'):
    ax.spines[sp].set_color('#8D99AE')
ax.tick_params(colors='#C9CCD8')
ax.set_xlabel('steps before reaching 1', color='#E6E6EE')
ax.set_ylabel(r'orbits of $1\leq n\leq 2000$, sorted by length', color='#E6E6EE')
ax.set_xticks([-175, -150, -125, -100, -75, -50, -25, 0])
ax.set_xticklabels(['175', '150', '125', '100', '75', '50', '25', '0'])
ax.set_yticks([])
ax.annotate('9232', (-34, 0.5), xytext=(-34, -30), ha='center', va='bottom', color='#FFF3C4', fontsize=7,
            annotation_clip=False, arrowprops=dict(arrowstyle='-', color='#FFF3C4', lw=0.5))
cax = fig.add_axes([0.08, 0.935, 0.89, 0.016])
cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=EMBER), cax=cax, orientation='horizontal')
cb.outline.set_visible(False)
cax.xaxis.set_ticks_position('top'); cax.xaxis.set_label_position('top')
cb.ax.tick_params(colors='#C9CCD8', labelsize=7)
cb.set_ticks([0, 4, 8, 12, 16, 20])
cb.set_label(r'$\log_2$ of the iterate', color='#E6E6EE', fontsize=8)
save(fig, 'art_g_heatmap', facecolor=NIGHT, bbox_inches=None)
