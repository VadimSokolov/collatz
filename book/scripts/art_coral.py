"""The Collatz coral: union of the orbits of 1..N drawn as turtle paths.

Each orbit is read backwards, from 1 out to its starting value, using the
standard map C(x) = x/2 (x even), 3x+1 (x odd).  Every value m != 1 has a
unique parent C(m); the segment from parent to m has unit length and heading
heading(C(m)) + ANG_EVEN if m is even, heading(C(m)) - ANG_ODD if m is odd.
Shared orbit tails coincide exactly, so the picture is the inverse Collatz
tree restricted to the orbits of 1..N, drawn with a fixed turning rule.
Each segment is coloured by the depth of m, i.e. its total stopping time.
"""
from art_common import *
from matplotlib.collections import LineCollection


def coral(N, ang_even=0.135, ang_odd=0.27, shrink=0.0, head0=np.pi / 2):
    pos = {1: (0.0, 0.0)}
    head = {1: head0}
    depth = {1: 0}
    segs, dep, child = [], [], []
    for n in range(2, N + 1):
        if n in pos:
            continue
        path = []
        x = n
        while x not in pos:
            path.append(x)
            x = 3 * x + 1 if x & 1 else x // 2
        for m in reversed(path):
            p = 3 * m + 1 if m & 1 else m // 2
            h = head[p] + (ang_even if m % 2 == 0 else -ang_odd)
            L = 1.0 if shrink == 0 else 1.0 / (1 + shrink * depth[p])
            x0, y0 = pos[p]
            q = (x0 + L * np.cos(h), y0 + L * np.sin(h))
            pos[m] = q; head[m] = h; depth[m] = depth[p] + 1
            segs.append(((x0, y0), q)); dep.append(depth[m]); child.append(m)
    return np.array(segs), np.array(dep), np.array(child), pos, depth


def draw_coral(ax, segs, dep, cmap, lw=(0.9, 0.12), alpha=0.55, vmax=None,
               gamma=1.0, crange=(0.98, 0.2)):
    vmax = vmax or dep.max()
    t = np.clip(dep / vmax, 0, 1) ** gamma
    order = np.argsort(dep)[::-1]          # shallow (trunk) drawn on top
    c0, c1 = crange
    cols = plt.get_cmap(cmap)(c0 + (c1 - c0) * t[order])
    cols[:, 3] = alpha
    widths = lw[0] * (1 - t[order]) ** 1.6 + lw[1]
    lc = LineCollection(segs[order], colors=cols, linewidths=widths,
                        capstyle='round')
    ax.add_collection(lc)
    ax.autoscale_view()
    ax.set_aspect('equal')
    ax.axis('off')
    return lc


if __name__ == '__main__':
    import time
    t = time.time()
    s, d, c, pos, depth = coral(20000)
    print(len(s), d.max(), time.time() - t)


def fit_frame(ax, segs, W, H, pad=0.04, shift=(0.0, 0.0)):
    """Make the axes fill a W x H inch figure and centre the drawing in it."""
    ax.set_position([0, 0, 1, 1])
    xs = segs[:, :, 0].ravel(); ys = segs[:, :, 1].ravel()
    x0, x1, y0, y1 = xs.min(), xs.max(), ys.min(), ys.max()
    cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    w, h = (x1 - x0) * (1 + 2 * pad), (y1 - y0) * (1 + 2 * pad)
    if w / h > W / H:
        h = w * H / W
    else:
        w = h * W / H
    cx -= shift[0] * w; cy -= shift[1] * h
    ax.set_xlim(cx - w / 2, cx + w / 2); ax.set_ylim(cy - h / 2, cy + h / 2)


from matplotlib.colors import LinearSegmentedColormap
EMBER = LinearSegmentedColormap.from_list(
    'ember', ['#1B1F5E', PAL['violet'], PAL['plum'], PAL['coral'], PAL['gold'], '#FFF3C4'])
LAGOON = LinearSegmentedColormap.from_list(
    'lagoon', ['#14123A', '#2C3E91', PAL['teal'], '#3FB8A9', '#A8E6CF', '#F7FFF0'])
