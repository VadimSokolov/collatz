"""Figures for Chapter 4 (Lyapunov functions and the limits of descent)."""
import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from collections import Counter
import sys as _s; _s.path.insert(0, BOOK + '/scripts'); from bookstyle import *
OUT = BOOK + "/figures/"
L2 = np.log(2)
def syr(x):
    m = 3*x+1; k = 0
    while m % 2 == 0: m //= 2; k += 1
    return m, k
def T(n): return n//2 if n % 2 == 0 else (3*n+1)//2
def T0(n): return n//2 if n % 2 == 0 else 3*n+1

cols = ["#c0392b", PAL['teal'], "#239b56", "#b9770e", "#7d3c98", "#566573"]

# ---------------- Figure A: stride increments coloured by k, and ascent runs
fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.9))
xs = np.arange(3, 20001, 2)
inc, ks = [], []
for x in xs:
    y, k = syr(int(x)); inc.append(np.log(y/x)); ks.append(k)
inc, ks = np.array(inc), np.array(ks)
xx = np.geomspace(3, 20000, 200)
for k in range(1, 7):
    sel = (ks == k) if k < 6 else (ks >= 6)
    lab = f"$k={k}$" if k < 6 else "$k\\geq 6$"
    ax[0].scatter(xs[sel], inc[sel], s=1.2, color=cols[k-1], label=lab, rasterized=True, lw=0)
    if k < 6:
        ax[0].plot(xx, np.log(3+1/xx) - k*L2, color=cols[k-1], lw=0.6)
# every k that occurs (up to k=14, at x=5461) gets its own thin curve, so that each point lies on the curve of its k
for k in range(6, ks.max() + 1):
    ax[0].plot(xx, np.log(3+1/xx) - k*L2, color=cols[5], lw=0.4)
ax[0].axhline(0, color="k", lw=0.5)
ax[0].axhline(np.log(3/4), color="k", lw=0.6, ls="--")
ax[0].text(3.5, np.log(3/4)+0.12, "mean $\\log(3/4)$", fontsize=7)
ax[0].set_xscale("log"); ax[0].set_xlabel("odd $x$"); ax[0].set_ylabel("$\\log S(x) - \\log x$")
# no point with k >= 9 has x < 341, so the lower left corner is free of data; an opaque legend hides only curves there
ax[0].legend(fontsize=6.5, markerscale=5, loc="lower left", ncol=2,
             frameon=True, facecolor="white", edgecolor="none", framealpha=1.0)
ax[0].set_ylim(-9.0, 0.75)
print("k range", ks.min(), ks.max(), "min increment", inc.min(), "at x =", xs[inc.argmin()])

# chained ascents from x = 2^j - 1
# grayscale-safe key: the four hues print at similar grays, so each orbit also gets its own dash pattern and marker (dash, marker, marker size)
ORB_STYLE = [("-", "o", 2.8), ((0, (6, 2.5)), "s", 2.5), ((0, (6, 2, 1, 2)), "^", 3.2), ((0, (1, 2)), "D", 2.2)]
for i, j in enumerate([6, 12, 18, 24]):
    x0 = 2**j - 1; x = x0; path = [0.0]
    for _ in range(int(2.2*j) + 12):
        if x == 1: break
        x, k = syr(x); path.append(np.log(x) - np.log(x0))
    ax[1].plot(range(len(path)), path, ls=ORB_STYLE[i][0], marker=ORB_STYLE[i][1], ms=ORB_STYLE[i][2], markevery=(i, 4), lw=1.0,
               color=cols[i], label=f"$x=2^{{{j}}}-1$")
jj = np.arange(0, 24)
ax[1].plot(jj, jj*np.log(1.5), "k--", lw=0.6, label="slope $\\log(3/2)$")
ax[1].axhline(0, color="gray", lw=0.4)
ax[1].set_xlabel("stride $i$"); ax[1].set_ylabel("$\\log S^{i}(x) - \\log x$")
ax[1].set_xlim(0, 60)
hs = [plt.Line2D([], [], color=cols[i], ls=ORB_STYLE[i][0], marker=ORB_STYLE[i][1], ms=ORB_STYLE[i][2], lw=1.0, label=f"$x=2^{{{j}}}-1$")
      for i, j in enumerate([6, 12, 18, 24])]
hs.append(plt.Line2D([], [], color="k", ls="--", lw=0.6, label="slope $\\log(3/2)$"))
ax[1].legend(handles=hs, frameon=False, fontsize=6.5, loc="upper right", handlelength=3.2)  # strides 40 to 60 stay below 5, so the upper right is free
plt.tight_layout()
plt.savefig(OUT + "ch4_stride_runs.pdf", bbox_inches="tight", dpi=300)

# ---------------- Figure B: clean redraw of the two step and phase portrait figure
fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.9))
odd = np.arange(3, 1001, 2)
two = np.array([np.log(T(T(int(x)))/x) for x in odd])
up = odd % 4 == 3
ax[0].scatter(np.log(odd[up]), two[up], s=3, color=cols[0], lw=0, label="$x\\equiv 3 \\ (\\mathrm{mod}\\ 4)$")
ax[0].scatter(np.log(odd[~up]), two[~up], s=3, color=cols[1], lw=0, label="$x\\equiv 1 \\ (\\mathrm{mod}\\ 4)$")
lx = np.linspace(np.log(3), np.log(1000), 200); xx = np.exp(lx)
ax[0].plot(lx, np.log(3/4 + 1/(4*xx)), "k--", lw=0.7, label="$\\log(3/4+1/(4x))$")
ax[0].plot(lx, np.log(9/4 + 5/(4*xx)), "k:", lw=0.8, label="$\\log(9/4+5/(4x))$")
ax[0].axhline(0, color="gray", lw=0.4)
ax[0].set_xlabel("$\\log x$"); ax[0].set_ylabel("$\\log T^{2}(x) - \\log x$")
ax[0].legend(fontsize=6.5, frameon=False, loc="center right")
cm = plt.cm.viridis
for n0 in [27, 97, 255, 703, 6171]:
    orb = [n0]
    while orb[-1] != 1: orb.append(T(orb[-1]))
    a = np.log(np.array(orb[:-1], float)); b = np.log(np.array(orb[1:], float))
    t = np.arange(len(a))/max(1, len(a)-1)
    for j in range(len(a)-1):
        ax[1].plot(a[j:j+2], b[j:j+2], color=cm(t[j]), lw=0.6, alpha=0.8)
# the accelerated orbits of 703 and 6171 peak at 125252 and 487700 (log 11.74 and 13.10), so the frame goes to 13.5
xr = np.linspace(0, 13.5, 50)
ax[1].plot(xr, xr - L2, color=cols[1], lw=1.0, label="even branch $y=x-\\log 2$")
ax[1].plot(xr, xr + np.log(1.5), color=cols[0], lw=1.0, label="odd branch $y=x+\\log(3/2)$")
ax[1].plot(xr, xr, "k:", lw=0.5)
ax[1].set_xlim(0, 13.5); ax[1].set_ylim(0, 13.5)
ax[1].set_xlabel("$\\log x_n$"); ax[1].set_ylabel("$\\log x_{n+1}$")
ax[1].legend(fontsize=6.5, frameon=False, loc="upper left")
plt.tight_layout()
plt.savefig(OUT + "ch4_phase.pdf", bbox_inches="tight", dpi=300)
