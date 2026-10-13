import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
# Figures for chapters 1 and 2 (shared book style)
#   ch1_orbit27.pdf, ch1_totalstop.pdf, ch2_residue_tree.pdf
import sys; sys.path.insert(0, BOOK + '/scripts'); from bookstyle import *
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

FIG = BOOK + "/figures/"
def T0(n): return n // 2 if n % 2 == 0 else 3 * n + 1
def T(n): return n // 2 if n % 2 == 0 else (3 * n + 1) // 2

# ---- ch1: orbit of 27 (odd iterates filled circles, even iterates open squares, so the parity survives in grayscale)
o = [27]
while o[-1] != 1: o.append(T0(o[-1]))
fig, ax = plt.subplots(figsize=(WIDTH, 2.5))
x = np.arange(len(o)); y = np.array(o)
ax.fill_between(x, 1, y, color=PAL['sky'], alpha=0.10, lw=0)
ax.plot(x, y, color=PAL['grey'], lw=0.8, zorder=1)
odd = y % 2 == 1
ax.scatter(x[odd], y[odd], s=16, color=PAL['coral'], zorder=3, label="odd iterate", lw=0)
ax.scatter(x[~odd], y[~odd], s=11, marker="s", facecolors="white", edgecolors=PAL['sky'], linewidths=0.9, zorder=2, label="even iterate")
ax.set_yscale("log")
k = int(np.argmax(y))
ax.annotate("peak 9232 at step %d" % k, (k, y[k]), xytext=(k - 62, y[k] * 1.25), fontsize=8,
            color=PAL['ink'], arrowprops=dict(arrowstyle="->", lw=0.7, color=PAL['ink']))
ax.set_xlabel("step $j$"); ax.set_ylabel("$T_0^{(j)}(27)$")
ax.set_xlim(-2, len(o) + 1)
ax.legend(loc="lower left")
fig.savefig(FIG + "ch1_orbit27.pdf"); plt.close(fig)

# ---- ch1: total stopping times for n <= 1e5 (red line = model prediction)
N = 10 ** 5
tst = np.zeros(N + 1, dtype=int)
for n in range(2, N + 1):
    m = n; c = 0
    while m >= n:
        m = T0(m); c += 1
    tst[n] = c + tst[m]
fig, ax = plt.subplots(figsize=(WIDTH, 2.7))
ns = np.arange(1, N + 1)
ax.scatter(ns, tst[1:], s=1.0, c=tst[1:], cmap=CMAP_SEQ, alpha=0.6, rasterized=True, linewidths=0)
ax.set_xscale("log")
xx = np.geomspace(2, N, 400)
ax.plot(xx, 3 * np.log(xx) / math.log(4 / 3), color=PAL['coral'], lw=1.6,
        label=r"$3\log n/\log(4/3)$")
ax.set_xlabel("$n$"); ax.set_ylabel("total stopping time under $T_0$")
ax.legend(loc="upper left")
fig.savefig(FIG + "ch1_totalstop.pdf"); plt.close(fig)
print("max tst", tst.max(), int(np.argmax(tst)))

# ---- ch2: stopping tree to depth 7 (red leaves stop, blue nodes continue)
def parvec(n, k):
    v = []
    for _ in range(k): v.append(n % 2); n = T(n)
    return v
def stopped(r, j):
    m = 0
    for i, b in enumerate(parvec(r, j), 1):
        m += b
        if 3 ** m < 2 ** i: return i
    return None
D = 7
nodes = {}; edges = []; pos = {}
layers = [[0]]
for j in range(1, D + 1):
    new = []
    for r in layers[-1]:
        for c in (r, r + 2 ** (j - 1)):
            new.append(c); edges.append(((j - 1, r), (j, c)))
    layers.append(new)
    nodes.update({(j, c): (stopped(c, j) == j) for c in new})
children = {}
for (a, b) in edges: children.setdefault(a, []).append(b)
def kids(nd):
    if nd != (0, 0) and nodes[nd]: return []
    return children.get(nd, [])
ycount = [0]
def place(nd):
    ks = kids(nd)
    if not ks or nd[0] == D:
        pos[nd] = ycount[0]; ycount[0] += 1; return pos[nd]
    ys = [place(k) for k in ks]; pos[nd] = sum(ys) / len(ys); return pos[nd]
place((0, 0))
fig, ax = plt.subplots(figsize=(WIDTH, 4.3))
ax.grid(False)
def draw(nd):
    for k in kids(nd):
        if nd[0] < D:
            ax.plot([nd[0], nd[0] + 0.5, k[0]], [pos[nd], pos[k], pos[k]], color=PAL['grey'],
                    lw=0.7, zorder=1)
            draw(k)
draw((0, 0))
for nd, yv in pos.items():
    j, r = nd
    if j == 0: lab = "all"; col = PAL['ink']
    else:
        lab = str(r)
        col = PAL['coral'] if nodes[nd] else PAL['teal']
    ax.text(j, yv, lab, ha="center", va="center", fontsize=6.5, color="white",
            bbox=dict(boxstyle="round,pad=0.25", fc=col, ec="none"), zorder=2)
ax.set_xticks(range(D + 1)); ax.set_xticklabels(["$2^{%d}$" % j for j in range(D + 1)])
ax.set_xlabel("modulus"); ax.set_yticks([]); ax.spines["left"].set_visible(False)
ax.invert_yaxis(); ax.set_xlim(-0.5, D + 0.5)
fig.savefig(FIG + "ch2_residue_tree.pdf"); plt.close(fig)
