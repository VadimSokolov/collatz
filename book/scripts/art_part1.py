"""Part I opener: the Collatz coral of the orbits of 1..20000.

Unit steps read from 1 outward; +0.135 rad at an even value, -0.27 rad at an
odd value (map 3x+1, x/2).  Colour and width encode depth in the inverse tree
(total stopping time): bright trunk near 1, dark blue tips at the longest
orbits.  4.5in x 3.5in, raster lines at 400 dpi.
"""
from art_coral import *

segs, dep, child, pos, depth = coral(20000, 0.135, 0.27, head0=2.15)
print('segments', len(segs), 'max depth', dep.max())
fig = plt.figure(figsize=(4.5, 3.5), facecolor=NIGHT)
ax = fig.add_subplot(); ax.set_facecolor(NIGHT)
lc = draw_coral(ax, segs, dep, LAGOON, lw=(1.1, 0.12), alpha=0.5,
                vmax=np.quantile(dep, 0.99), gamma=0.8, crange=(0.99, 0.22))
lc.set_rasterized(True)
fit_frame(ax, segs, 4.5, 3.5, pad=0.03)
save(fig, 'art_part1', facecolor=NIGHT, bbox_inches=None)
