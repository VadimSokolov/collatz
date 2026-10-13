"""Figure for Chapter 8: density of the Erdos 125 sumset A+B in [0,N], N up to 4^13
(data from ch08_erdos125.py), in the shared book style."""
import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
import sys; sys.path.insert(0, BOOK + '/scripts'); from bookstyle import *
import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
g, d = np.load(BOOK + "/scripts/ch08_erdos125_curve.npy")
fig, ax = plt.subplots(figsize=(WIDTH, 2.9))
for m in range(7, 17):
    ax.axvline(3**m, c=PAL['violet'], lw=0.6, alpha=0.35, zorder=0)
for k in range(5, 14):
    ax.axvline(4**k, c=PAL['green'], lw=0.6, ls=':', alpha=0.6, zorder=0)
ax.fill_between(g, 1015/1458, d, color=PAL['sky'], alpha=0.18, lw=0)
ax.plot(g, d, c=PAL['teal'], lw=1.0, label=r"density of $A+B$ in $[0,N]$")
ax.axhline(1015/1458, c=PAL['coral'], ls="--", lw=1.1, label="Hasler and Melfi bound on the lower density")
i = int(np.argmin(np.where(g >= 1e4, d, 9)))
ax.plot(g[i], d[i], "o", ms=4, c=PAL['gold'], mec=PAL['ink'], mew=0.4, label=r"minimum, $N=3^{15}-1$")
ax.plot([], [], c=PAL['violet'], lw=0.6, alpha=0.6, label=r"$N=3^m$")
ax.plot([], [], c=PAL['green'], lw=0.6, ls=':', label=r"$N=4^k$")
ax.set_xscale("log"); ax.set_xlabel(r"$N$"); ax.set_ylabel("density")
ax.set_ylim(0.5, 1.0); ax.set_xlim(1e3, 4**13)
ax.legend(fontsize=6.3, loc="lower right", ncol=2)
fig.savefig(BOOK + "/figures/ch8_erdos125.pdf"); plt.close(fig)
