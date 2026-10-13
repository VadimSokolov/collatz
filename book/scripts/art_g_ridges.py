"""Gallery: ridgelines of the stopping time record orbits.

The records for total stopping time sigma_inf(n) (n with sigma_inf(n) larger
than for every smaller n) are computed for n <= 10^6.  For each record holder
n >= 27, one ridge shows log2 T_0^(j)(n) against the step j, from j = 0 to the
arrival at 1.  Ridges are stacked with the smallest record at the top and
filled to hide the ridges behind, in the manner of a joy plot (heights 0.42 units per unit of log2,
spacing 4 units); colour runs
with the record's rank.  Gold dots mark each orbit's maximum.
"""
from art_common import *
from art_coral import EMBER

N = 10 ** 6
st, mx = total_stopping(N)
recs = []; best = -1
for n in range(1, N + 1):
    if st[n] > best:
        best = st[n]; recs.append(n)
print('number of records up to 1e6:', len(recs))
print(recs)
R = [n for n in recs if n >= 27]
print('ridges:', len(R), 'from', R[0], 'to', R[-1], 'steps', st[R[0]], 'to', st[R[-1]])

fig = plt.figure(figsize=(4.5, 6.2), facecolor=NIGHT)
ax = fig.add_axes([0.04, 0.06, 0.92, 0.92]); ax.set_facecolor(NIGHT); ax.axis('off')
gap = 4.0
hs = 0.42   # vertical scale: one unit of log2 is 0.42 of a unit, ridge spacing 4 units
nR = len(R)
maxlog = 0
for i, n in enumerate(R):
    o = np.log2(np.array(collatz_orbit(n), float))
    maxlog = max(maxlog, o.max())
    o = hs * o
    base = (nR - 1 - i) * gap
    j = np.arange(len(o))
    col = EMBER(0.25 + 0.75 * i / (nR - 1))
    z = 2 * i
    ax.fill_between(j, base, base + o, color=NIGHT, lw=0, zorder=z)
    ax.plot(j, base + o, color=col, lw=0.8, zorder=z + 1)
    ax.plot([0, len(o) - 1], [base, base], color=col, lw=0.3, alpha=0.5, zorder=z + 1)
    k = np.argmax(o)
    ax.scatter([k], [base + o[k]], s=5, color='#FFF3C4', lw=0, zorder=z + 1)
    ax.text(-5, base + o[0], f'{n}', ha='right', va='center', fontsize=5.5, color=col, zorder=z + 1)
print('largest log2 of any iterate', maxlog)
ax.set_xlim(-40, 530)
ax.set_ylim(-3, (nR - 1) * gap + hs * 34)
for t in range(0, 501, 100):
    ax.text(t, -2, f'{t}', ha='center', va='top', fontsize=6.5, color='#C9CCD8')
ax.text(250, -7.5, r'step $j$', ha='center', va='top', fontsize=8, color='#E6E6EE')
save(fig, 'art_g_ridges', facecolor=NIGHT, bbox_inches=None)
