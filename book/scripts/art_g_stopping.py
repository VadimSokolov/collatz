"""Gallery: the stopping time landscape for 2 <= n <= 10^6.

Every n is one point (n, sigma_inf(n)), where sigma_inf(n) is the number of
steps of T_0 (x -> 3x+1 for odd x, x/2 for even x) needed to reach 1.  The
colour of the point is the excursion ratio log2(M(n)/n), M(n) the largest
value on the orbit.  Points are drawn in increasing order of excursion so the
rare high flyers sit on top.  Raster layer at 400 dpi.
"""
from art_common import *
from art_coral import EMBER

N = 10**6
st, mx = total_stopping(N)
n = np.arange(2, N + 1)
s = st[2:]; ex = np.log2(mx[2:] / n)
i = np.argmax(s)
print('max sigma_inf', s[i], 'at n =', n[i])
j = np.argmax(ex)
print('max log2 excursion', ex[j], 'at n =', n[j], 'M =', mx[2:][j])
print('fraction with ex < 1:', np.mean(ex < 1), ' median sigma', np.median(s))
print('fraction of n with M(n) = n (n power of 2 times ... ):', np.mean(ex == 0))

fig = plt.figure(figsize=(4.5, 3.3), facecolor=NIGHT)
ax = fig.add_axes([0.12, 0.13, 0.74, 0.83]); ax.set_facecolor(NIGHT)
o = np.argsort(ex, kind='stable')
CM = mpl.colors.LinearSegmentedColormap.from_list('emb2', EMBER(np.linspace(0.12, 1, 256)))
norm = mpl.colors.PowerNorm(0.6, 0, 17)
sc = ax.scatter(n[o], s[o], c=ex[o], cmap=CM, norm=norm, s=0.35, lw=0, alpha=0.7)
sc.set_rasterized(True)
ax.set_xlim(0, N); ax.set_ylim(0, 545)
ax.grid(False)
for sp in ('left', 'bottom'):
    ax.spines[sp].set_color('#8D99AE')
ax.tick_params(colors='#C9CCD8')
ax.set_xticks([0, 2.5e5, 5e5, 7.5e5, 1e6])
ax.set_xticklabels(['0', '250000', '500000', '750000', '1000000'])
ax.set_xlabel(r'starting value $n$', color='#E6E6EE')
ax.set_ylabel(r'total stopping time $\sigma_\infty(n)$', color='#E6E6EE')
ax.annotate(f'$n = {n[i]}$, ${s[i]}$ steps', (n[i], s[i]), xytext=(-60, 0), ha='right',
            textcoords='offset points', color='#FFF3C4', fontsize=7, va='center',
            arrowprops=dict(arrowstyle='-', color='#FFF3C4', lw=0.5, shrinkB=2))
ax.scatter([n[i]], [s[i]], s=18, facecolor='none', edgecolor='#FFF3C4', lw=0.7)
cax = fig.add_axes([0.89, 0.13, 0.018, 0.83])
cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=CM), cax=cax)
cb.outline.set_visible(False)
cb.ax.tick_params(colors='#C9CCD8', labelsize=7)
cb.set_ticks([0, 1, 2, 4, 6, 8, 11, 14, 17])
cb.set_label(r'excursion $\log_2\,(M(n)/n)$', color='#E6E6EE')
save(fig, 'art_g_stopping', facecolor=NIGHT, bbox_inches=None)
