"""Gallery: the inverse tree of 1 under T_0, to depth 20, as a radial tree.

Node m at depth d is a number with sigma_inf(m) = d.  Its children are 2m
(always) and (m-1)/3 when m = 4 mod 6 and m > 4 (the exception removes the
trivial cycle 1, 4, 2).  Radius = depth.  Every node below depth 20 has the
child 2m, so the leaves are exactly the numbers with sigma_inf = 20, spaced
evenly round the rim; each inner node sits at the mean angle of its leaves.
Edges are drawn as polar elbows (arc at the parent radius, then a radial
step).  Gold edges and rings mark the odd branches m -> (m-1)/3.
"""
from art_common import *
from art_coral import EMBER

D = 20
children = {}
depth = {1: 0}
lev = [1]; order = [1]
for d in range(1, D + 1):
    nx = []
    for m in lev:
        ch = [2 * m]
        if m % 6 == 4 and m > 4:
            ch.append((m - 1) // 3)
        children[m] = ch
        for c in ch:
            depth[c] = d
        nx += ch
    lev = nx
leaves = lev
counts = {d: sum(1 for v in depth.values() if v == d) for d in range(D + 1)}
print('nodes', len(depth), 'leaves', len(leaves), 'level sizes', [counts[d] for d in range(D + 1)])
odd = [m for m in depth if m % 2 == 1]
print('odd nodes', len(odd), 'largest leaf', max(leaves), 'smallest leaf', min(leaves))

# angular layout by depth first leaf order (odd branch first, then doubling)
ang = {}
k = [0]
def place(m):
    if depth[m] == D:
        ang[m] = 2 * np.pi * (k[0] + 0.5) / len(leaves); k[0] += 1
        return ang[m]
    ch = sorted(children[m], key=lambda c: c % 2 == 0)
    a = [place(c) for c in ch]
    ang[m] = np.mean(a)
    return ang[m]
place(1)
rot = np.pi / 2 - ang[1]
for m in ang:
    ang[m] += rot

fig = plt.figure(figsize=(4.5, 4.5), facecolor=NIGHT)
ax = fig.add_axes([0, 0, 1, 1]); ax.set_facecolor(NIGHT); ax.axis('off')
ax.set_aspect('equal')
R = D + 5.2
ax.set_xlim(-R, R); ax.set_ylim(-R, R)
for r in range(1, D + 1):
    t = np.linspace(0, 2 * np.pi, 400)
    ax.plot(r * np.cos(t), r * np.sin(t), color='#2A2950', lw=0.3, zorder=0)
cm = lambda d: EMBER(0.72 - 0.52 * d / D)
for m, ch in children.items():
    r0 = depth[m]; a0 = ang[m]
    for c in ch:
        a1 = ang[c]
        col = PAL['gold'] if c % 2 else cm(depth[c])
        lw = 1.1 if c % 2 else 0.8
        t = np.linspace(a0, a1, max(2, int(abs(a1 - a0) * 60)))
        xs = np.r_[r0 * np.cos(t), (r0 + 1) * np.cos(a1)]
        ys = np.r_[r0 * np.sin(t), (r0 + 1) * np.sin(a1)]
        ax.plot(xs, ys, color=col, lw=lw, alpha=0.95, solid_capstyle='round', zorder=2)
for m in depth:
    r = depth[m]; a = ang[m]
    x, y = r * np.cos(a), r * np.sin(a)
    if m % 2:
        ax.scatter(x, y, s=9, facecolor=NIGHT, edgecolor=PAL['gold'], lw=0.8, zorder=4)
    else:
        ax.scatter(x, y, s=4, color=cm(r), lw=0, zorder=3)
# labels: leaves radially outside the rim
for m in leaves:
    a = ang[m]; deg = np.degrees(a) % 360
    flip = 90 < deg < 270
    ax.text((D + 0.6) * np.cos(a), (D + 0.6) * np.sin(a), f'{m}', rotation=deg + (180 if flip else 0),
            rotation_mode='anchor', ha='right' if flip else 'left', va='center', fontsize=4.5,
            color='#E8D9F5' if m % 2 == 0 else '#FFE3A3')
# labels of inner nodes on the spine of small depths
for m in depth:
    d = depth[m]
    if 0 < d <= 9 or m in (1,):
        a = ang[m]
        tx, ty = -np.sin(a), np.cos(a)
        ax.text(d * np.cos(a) + 0.62 * tx, d * np.sin(a) + 0.62 * ty, f'{m}', fontsize=4.8,
                color='#FFF3C4', ha='center', va='center', zorder=5)
save(fig, 'art_g_tree', facecolor=NIGHT, bbox_inches=None)
