"""Figures for Chapter 8: records of r/(3 Lambda_r) (ch8_records.pdf) and one sided
distances (ch8_onesided.pdf), in the shared book style."""
import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
import sys; sys.path.insert(0, BOOK + '/scripts'); from bookstyle import *
import mpmath as mp, numpy as np, math
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

mp.mp.dps = 120
TH = mp.log(3)/mp.log(2); L2 = mp.log(2); L3 = mp.log(3)
D = 10**100; T = int(mp.floor(TH*D))

def cf(x, n):
    a = []; y = mp.mpf(x)
    for _ in range(n):
        ai = int(mp.floor(y)); a.append(ai); y = 1/(y-ai)
    return a
A = cf(TH, 40); q = [0, 1]
for ai in A[1:32]:
    q.append(ai*q[-1]+q[-2])
Q = lambda k: q[k+1]
def fval(r):
    s = int(mp.ceil(r*TH)); return float(r/(3*(s*L2-r*L3)))

# ---------- records of F
N = 10**5
r = np.arange(1, N+1)
dd = np.array([((-(k*T)) % D)/D for k in range(1, N+1)], dtype=float)
fall = r/(3*dd*math.log(2))
upper = []
for k in range(0, 30, 2):
    for j in range(1, A[k+1]+1):
        upper.append(Q(k-1)+j*Q(k))
upper = sorted(set(u for u in upper if u <= 10**12))
fu = [fval(u) for u in upper]
fig, ax = plt.subplots(figsize=(WIDTH, 3.2))
ax.scatter(r, fall, s=0.5, c=np.log10(dd), cmap=CMAP_HOT, vmin=-5, vmax=0.3, rasterized=True, lw=0, alpha=0.8)
ax.plot([], [], "o", c=PAL['plum'], ms=2.5, label=r"all $r\leq 10^5$ (shade: $\log\delta(r)$)")
ax.step(upper, fu, where="post", c=PAL['teal'], lw=1.0)
ax.plot(upper, fu, "o", ms=3.0, c=PAL['teal'], mec=PAL['ink'], mew=0.3, label="one sided ladder rungs")
cu = [Q(k) for k in range(1, 25, 2) if Q(k) <= 10**12]
ax.plot(cu, [fval(c) for c in cu], "s", mfc="none", mec=PAL['coral'], mew=1.0, ms=5.5, label="convergents from above")
cl = [Q(k) for k in range(2, 25, 2) if Q(k) <= 10**12]
ax.plot(cl, [fval(c) for c in cl], "v", c=PAL['green'], mec=PAL['ink'], mew=0.3, ms=4.5, label="convergents from below")
ax.axhline(2.0**68, c=PAL['gold'], ls="--", lw=1.2); ax.text(1.2, 2.0**68*4, r"$2^{68}$", fontsize=8)
ax.annotate(r"$r=72\,057\,431\,991$", xy=(72057431991, fval(72057431991)), xytext=(1e6, 1e22),
            fontsize=7, arrowprops=dict(arrowstyle="->", lw=0.6, color=PAL['ink']))
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel(r"number of odd elements $r$"); ax.set_ylabel(r"$F(r)=r/(3\Lambda_r)$")
ax.set_ylim(0.3, 1e24); ax.legend(fontsize=6.5, loc="lower right", markerscale=1.0)
fig.savefig(BOOK + "/figures/ch8_records.pdf"); plt.close(fig)

# ---------- one sided distance
N2 = 700
r2 = np.arange(1, N2+1); d2 = dd[:N2]
two = np.minimum(d2, 1-d2)
mins = []; m = 2
for k, v in zip(r2, d2):
    if v < m: m = v; mins.append(k)
fig, ax = plt.subplots(figsize=(WIDTH, 2.9))
ax.scatter(r2, d2, s=3, c=r2 % 53, cmap=CMAP_CYC, lw=0, alpha=0.75,
           label=r"one sided distance $\lceil r\theta\rceil-r\theta$")
ax.step(list(mins)+[N2], list(d2[np.array(mins)-1])+[d2[mins[-1]-1]], where="post", c=PAL['teal'], lw=1.0)
ax.plot(mins, d2[np.array(mins)-1], "o", ms=4.5, mfc="white", mec=PAL['teal'], mew=1.2,
        label="one sided record minima")
for c in [12, 53, 665]:
    ax.plot(c, d2[c-1], "v", c=PAL['coral'], ms=5, mec=PAL['ink'], mew=0.3)
    ax.plot(c, two[c-1], "v", mfc="none", mec=PAL['coral'], mew=1.0, ms=5)
    ax.plot([c, c], [two[c-1], d2[c-1]], c=PAL['coral'], lw=0.6, ls=':')
ax.plot([], [], "v", c=PAL['coral'], ms=5, label="convergents from below, one sided")
ax.plot([], [], "v", mfc="none", mec=PAL['coral'], ms=5, label="convergents from below, two sided")
ax.set_yscale("log"); ax.set_xlabel(r"$r$"); ax.set_ylabel("distance")
ax.set_ylim(3e-5, 3); ax.legend(fontsize=6.5, loc="lower left")
fig.savefig(BOOK + "/figures/ch8_onesided.pdf"); plt.close(fig)
print("mins", mins)
