"""Front cover: a clear Collatz coral of the orbits of 1..N (no title text).

Turning rule: unit steps read from 1 outward, +ANG_EVEN rad at an even value,
-ANG_ODD rad at an odd value (standard map 3x+1, x/2).  Colour encodes depth
in the inverse tree (total stopping time); width decreases with depth.
Fewer orbits and bolder strokes than the plates, so the structure reads from
across a room.  Raster layer at 400 dpi inside a 7.25in x 10.5in PDF: the front panel
of the 7 x 10 in paperback cover with its bleed.
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from art_coral import *
import sys
N = int(sys.argv[1]) if len(sys.argv) > 1 else 2500
AE, AO = (float(sys.argv[2]), float(sys.argv[3])) if len(sys.argv) > 3 else (0.17, 0.30)
H0 = float(sys.argv[4]) if len(sys.argv) > 4 else 0.45
segs, dep, child, pos, depth = coral(N, AE, AO, head0=H0)
print('N', N, 'segments', len(segs), 'max depth', dep.max())
PW, PH = 7.25, 10.5                          # page size in inches
fig = plt.figure(figsize=(PW, PH), facecolor=NIGHT)
ax = fig.add_subplot(); ax.set_facecolor(NIGHT)
lc = draw_coral(ax, segs, dep, EMBER, lw=(3.8, 0.55), alpha=0.92,
                vmax=np.quantile(dep, 0.985), gamma=0.8, crange=(0.99, 0.12))
lc.set_rasterized(True)
fit_frame(ax, segs, PW, PH, pad=0.05)
# leave the upper third for the title: expand the view upward, keeping the page aspect
x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim(); h = y1 - y0; cx = 0.5 * (x0 + x1)
top = 0.30                                   # fraction of the page reserved above the drawing
H = h / (1 - top); W = H * PW / PH
ax.set_ylim(y0 - 0.02 * h, y0 - 0.02 * h + H); ax.set_xlim(cx - W / 2, cx + W / 2)
save(fig, 'art_cover', facecolor=NIGHT, bbox_inches=None)
