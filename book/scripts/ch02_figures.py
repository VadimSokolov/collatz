"""New and regenerated figures for Chapter 2 (shared book style).
   ch2_syracuse.pdf    Geometric law of v_2(3n+1) and the Syracuse increment
   ch2_sunburst.pdf    Terras parity tree of residues mod 2^k as a sunburst
   ch2_stopheat.pdf    stopping time sigma(n) arranged by residue class mod 64
   ch2_adic_shift.pdf  the map T in 2 adic coordinates and its conjugate, the shift
"""
import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
import sys; sys.path.insert(0, BOOK + '/scripts'); from bookstyle import *
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
import matplotlib.patheffects as pe
from matplotlib.colors import Normalize, LogNorm
from collections import Counter

FIG = BOOK + "/figures/"
def T(n): return n // 2 if n % 2 == 0 else (3 * n + 1) // 2
def v2(n):
    c = 0
    while n % 2 == 0: n //= 2; c += 1
    return c

# ------------------------------------------------------------------ ch2_syracuse
N = 2 ** 20
ks = np.array([v2(3 * n + 1) for n in range(1, N, 2)])
cnt = Counter(ks.tolist()); K = 10
xs = np.arange(1, K + 1)
emp = np.array([cnt[k] for k in xs]) / len(ks)
fig, axs = plt.subplots(1, 2, figsize=(WIDTH, 2.2))
axs[0].bar(xs - 0.2, emp, width=0.4, color=PAL['teal'], label="odd $n<2^{20}$")
axs[0].bar(xs + 0.2, 2.0 ** (-xs), width=0.4, color=PAL['gold'], label="$2^{-k}$")
axs[0].set_yscale("log"); axs[0].set_xlabel("$k=v_2(3n+1)$"); axs[0].set_ylabel("proportion")
axs[0].set_xticks(xs); axs[0].set_ylim(5e-4, 3); axs[0].legend(fontsize=7, ncol=2, loc="upper center", columnspacing=0.8)
vals = np.array([math.log(3) - k * math.log(2) for k in xs])
axs[1].bar(vals, emp, width=0.13, color=[PAL['coral'] if v > 0 else PAL['teal'] for v in vals])
axs[1].axvline(math.log(3 / 4), color=PAL['ink'], lw=0.9, ls="--")
axs[1].text(math.log(3 / 4) - 0.12, 0.42, r"mean $\log(3/4)$", fontsize=7, ha="right")
axs[1].text(math.log(1.5) - 0.12, 0.47, r"$\log(3/2)$", fontsize=7, color=PAL['coral'], ha="right")
axs[1].set_xlabel(r"$\log 3-k\log 2$"); axs[1].set_ylabel("proportion")
fig.tight_layout(); fig.savefig(FIG + "ch2_syracuse.pdf"); plt.close(fig)

# ------------------------------------------------------------------ ch2_sunburst
# Angular width of a node is proportional to the number of leaves below it, where the
# leaves are the classes that stop (tau = k) and the classes still alive at depth DEPTH.
DEPTH = 9
def residue_of_word(w):
    k = len(w); c = 0; m = 0
    for i, b in enumerate(w):
        if b: c = 3 * c + 2 ** i; m += 1
    return (-c * pow(3, -m, 2 ** k)) % 2 ** k
def height(w):
    return sum(w) * math.log(3) - len(w) * math.log(2)
def leaves(w):
    if len(w) > 0 and height(w) < 0: return 1
    if len(w) == DEPTH: return 1
    return leaves(w + (0,)) + leaves(w + (1,))
fig = plt.figure(figsize=(WIDTH, 3.9))
ax = fig.add_axes([0.0, 0.0, 0.80, 1.0]); ax.set_aspect("equal"); ax.axis("off")
cax = fig.add_axes([0.83, 0.2, 0.025, 0.6])
norm = Normalize(vmin=0.0, vmax=3.2)
cmap = plt.get_cmap(CMAP_SEQ)
r0, dr = 0.9, 0.42
ax.add_patch(plt.Circle((0, 0), r0 - 0.03, color=PAL['ink']))
ax.text(0, 0, "all", color="white", ha="center", va="center", fontsize=8)
TOT = leaves(())
def draw(w, a0, a1):
    k = len(w)
    if k > 0:
        h = height(w); rin = r0 + (k - 1) * dr
        if h < 0:  # stopped class: white wedge with coral hatch and outline, so that it prints without relying on hue
            ax.add_patch(Wedge((0, 0), rin + dr, 90 - a1, 90 - a0, width=dr, fc="white", ec=PAL['coral'],
                               hatch="////", lw=0.6))
        else:
            ax.add_patch(Wedge((0, 0), rin + dr, 90 - a1, 90 - a0, width=dr, fc=cmap(norm(h)), ec="white",
                               lw=0.5 if k > 5 else 0.8))
        mid = math.radians(90 - (a0 + a1) / 2)
        res = residue_of_word(w)
        if k <= 3 or h < 0:
            rr = rin + dr / 2
            rot = (90 - (a0 + a1) / 2)
            if k > 3:
                rot = rot if -90 <= rot <= 90 else rot - 180
            else:
                rot = 0
            ax.text(rr * math.cos(mid), rr * math.sin(mid), str(res), ha="center", va="center",
                    fontsize=5.6 if k > 3 else 7, color=PAL['ink'] if h < 0 else "white", rotation=rot,
                    path_effects=[pe.withStroke(linewidth=1.4, foreground="white")] if h < 0 else [])
        if k == DEPTH and h >= 0:
            rr = rin + dr + 0.28
            rot = 90 - (a0 + a1) / 2
            rot2 = rot if -90 <= rot <= 90 else rot - 180
            ax.text(rr * math.cos(mid), rr * math.sin(mid), str(res), ha="center", va="center",
                    fontsize=5.2, color=PAL['ink'], rotation=rot2)
        if h < 0 or k == DEPTH: return
    n0 = leaves(w + (0,)); n1 = leaves(w + (1,))
    am = a0 + (a1 - a0) * n0 / (n0 + n1)
    draw(w + (0,), a0, am); draw(w + (1,), am, a1)
draw((), 0.0, 360.0)
R = r0 + DEPTH * dr + 0.55
ax.set_xlim(-R, R); ax.set_ylim(-R, R)
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
cb = fig.colorbar(sm, cax=cax)
cb.set_label(r"height $m_k\log 3-k\log 2$ of a live class", fontsize=7)
cb.ax.tick_params(labelsize=6.5)
# key for stopped classes, no colour named: a hatched swatch, then the old text minus its first word (same right edge, so the page box is unchanged)
fig.add_artist(plt.Rectangle((0.800, 0.877), 0.049, 0.024, transform=fig.transFigure, fc="white", ec=PAL['coral'], hatch="////", lw=0.6))
fig.text(0.8694, 0.88, "stops at ring $k$", fontsize=7, color=PAL['ink'])
with plt.rc_context({"hatch.linewidth": 0.45}):
    fig.savefig(FIG + "ch2_sunburst.pdf")
plt.close(fig)
print("sunburst leaves", TOT)

# ------------------------------------------------------------------ ch2_stopheat
COLS, ROWS = 64, 128
sig = np.full(COLS * ROWS, np.nan)
for n in range(2, COLS * ROWS):
    x = n; c = 0
    while True:
        x = T(x); c += 1
        if x < n: break
    sig[n] = c
grid = sig.reshape(ROWS, COLS)
fig, ax = plt.subplots(figsize=(WIDTH, 3.3))
ax.grid(False)
cm = plt.get_cmap(CMAP_HOT).copy(); cm.set_bad(PAL['sand'])
im = ax.imshow(grid, aspect="auto", cmap=cm, norm=LogNorm(vmin=1, vmax=np.nanmax(grid)),
               interpolation="nearest", origin="upper")
ax.set_xlabel(r"column $n\ \mathrm{mod}\ 64$"); ax.set_ylabel(r"row $\lfloor n/64\rfloor$")
ax.set_xticks([0, 16, 32, 48, 63])
top = ax.secondary_xaxis("top")
top.set_xticks([7, 15, 27, 31, 39, 47, 59, 63]); top.tick_params(labelsize=6, length=2)
cb = fig.colorbar(im, ax=ax, pad=0.02)
cb.set_label(r"stopping time $\sigma(n)$")
cb.set_ticks([1, 2, 4, 8, 16, 32, 64]); cb.set_ticklabels(["1", "2", "4", "8", "16", "32", "64"])
fig.savefig(FIG + "ch2_stopheat.pdf"); plt.close(fig)
print("max sigma below 2^13:", np.nanmax(grid), int(np.nanargmax(sig)))
# columns that are constant (decided by n mod 64) versus columns that vary
const = [c for c in range(COLS) if np.nanmin(grid[:, c]) == np.nanmax(grid[:, c])]
print("constant columns", len(const), "varying columns", sorted(set(range(COLS)) - set(const)))

# ------------------------------------------------------------------ ch2_adic_shift
B = 16
M = 2 ** B
def iota(x):
    """2 adic coordinate in [0,1): x_0/2 + x_1/4 + ... (first B binary digits)."""
    x %= M; s = 0.0
    for i in range(B):
        s += ((x >> i) & 1) * 2.0 ** (-i - 1)
    return s
def parity_coord(x):
    s = 0.0
    for i in range(B):
        s += (x & 1) * 2.0 ** (-i - 1); x = T(x)
    return s
xs = np.arange(M)
X = np.array([iota(int(x)) for x in xs]); Y = np.array([iota(T(int(x))) for x in xs])
Xq = np.array([parity_coord(int(x)) for x in xs]); Yq = np.array([parity_coord(T(int(x))) for x in xs])
par = xs % 2 == 1
fig, axs = plt.subplots(1, 2, figsize=(WIDTH, 2.35))
for a, (u, v), ttl in zip(axs, [(X, Y), (Xq, Yq)],
                          [r"$T$ in $2$ adic coordinates", r"conjugate $Q\circ T\circ Q^{-1}=\mathrm{sh}$"]):
    a.scatter(u[~par], v[~par], s=0.6, color=PAL['teal'], lw=0, rasterized=True, label="even $x$")
    a.scatter(u[par], v[par], s=0.6, color=PAL['coral'], lw=0, rasterized=True, label="odd $x$")
    a.set_xlim(0, 1); a.set_ylim(0, 1); a.set_aspect("equal")
    a.set_title(ttl, fontsize=8)
    a.set_xticks([0, 0.5, 1]); a.set_yticks([0, 0.5, 1])
axs[0].set_xlabel(r"$\iota(x)$"); axs[0].set_ylabel(r"$\iota(T(x))$")
axs[1].set_xlabel(r"$\iota(Q(x))$"); axs[1].set_ylabel(r"$\iota(Q(T(x)))$")
axs[0].legend(loc="upper left", markerscale=8, fontsize=6.5, handletextpad=0.2)
fig.tight_layout(); fig.savefig(FIG + "ch2_adic_shift.pdf"); plt.close(fig)
print("done")
