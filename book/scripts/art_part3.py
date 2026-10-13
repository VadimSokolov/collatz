"""Part III opener: the sunflower of r log2 3 modulo 1.

For r = 1, ..., 3000 a point is placed at angle 2 pi {r theta} (theta = log2 3)
and radius sqrt(r) (Vogel's phyllotaxis layout, which spreads the points
uniformly over the disc).  The colour of each point is the one sided distance
delta(r) = ceil(r theta) - r theta, a circular coordinate, shown with a cyclic
colormap.  The spiral arms (parastichies) are residue classes of r modulo
continued fraction denominators of theta: nearest neighbours differ by 12, 41
or 53 near the centre and by 53 or 106 in the outer ring (53 theta = 84.0030).  Gold rings: the
delta records (r with delta(r) smaller than for every smaller r), which are the
rungs of the semiconvergent ladder; they approach the ray at angle 0 from the
clockwise side.  Computed in exact integer arithmetic via a 60 digit theta.
"""
from art_common import *
from mpmath import mp, mpf, log as mlog, floor as mfloor
from matplotlib.colors import LinearSegmentedColormap

mp.dps = 60
R = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
th = mlog(3) / mlog(2)
# fractional parts with 60 digit accuracy: frac(r th) = (r * TH mod 10^50) / 10^50
SCALE = 10**50
TH = int(mfloor(th * SCALE))
frac = np.array([((r * TH) % SCALE) / SCALE for r in range(1, R + 1)])
delta = 1 - frac                              # ceil(r th) - r th (never zero)
r = np.arange(1, R + 1)
recs = []; best = 2.0
for i in range(R):
    if delta[i] < best:
        best = delta[i]; recs.append(i + 1)
print('delta records up to', R, ':', recs)

ang = 2 * np.pi * frac
rad = np.sqrt(r)
X, Y = rad * np.cos(ang), rad * np.sin(ang)

SUN = LinearSegmentedColormap.from_list(
    'sun', [PAL['teal'], '#3FB8A9', PAL['gold'], PAL['coral'], PAL['plum'], PAL['violet'],
            '#2C3E91', PAL['teal']])
W, H = 4.5, 3.5
fig = plt.figure(figsize=(W, H), facecolor=NIGHT)
ax = fig.add_axes([0, 0, 1, 1]); ax.set_facecolor(NIGHT); ax.axis('off')
Rm = np.sqrt(R) * 1.02
ax.set_ylim(-Rm, Rm); ax.set_xlim(-Rm * W / H, Rm * W / H)
ax.set_aspect('equal')
sz = 9000 / R * 0.9 * (0.35 + 0.65 * rad / rad.max())
sc = ax.scatter(X, Y, c=delta, cmap=SUN, s=sz, lw=0, alpha=0.9, vmin=0, vmax=1)
sc.set_rasterized(True)
ax.plot([0, Rm], [0, 0], color='white', lw=0.5, alpha=0.5)
rr = np.array(recs)
ax.plot(X[rr - 1], Y[rr - 1], color='#FFF3C4', lw=0.6, alpha=0.7)
ax.scatter(X[rr - 1], Y[rr - 1], s=16, facecolor='none', edgecolor=PAL['gold'], lw=0.9, zorder=5)
ax.scatter(X[rr - 1], Y[rr - 1], s=3, color='#FFF3C4', lw=0, zorder=6)
for q in recs:
    if q >= 306:
        ax.annotate(f'{q}', (X[q - 1], Y[q - 1]), xytext=(4, -7), textcoords='offset points',
                    color='#FFF3C4', fontsize=6.5, fontstyle='italic')
save(fig, 'art_part3', facecolor=NIGHT, bbox_inches=None)
