"""Figure 20.1: (a) the odd elements of the Christoffel rational cycle at r = 15601,
as log2(x_t / M) against t; (b) the constant of the minimum inequality along the ladder."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sharp import christoffel_cycle
from thresholds import rows
import mpmath as mp
from ladder import TH
plt.rcParams.update({"font.family": "serif", "font.size": 9})
fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.6))
# panel (a): recompute x_t for r = 15601
r, s = 15601, 24727
t = np.arange(r); K = (t*s)//r
W = np.array([float(int(k) - int(x)*TH) for k, x in zip(K, t)])
a = np.exp2(W); d = s - r*TH
x0 = a.sum()/(3*float(mp.expm1(d*mp.log(2))))
x = (x0 + np.concatenate(([0.0], np.cumsum(a)[:-1]))/3)/a
M = x.min()
ax[0].plot(t[:120], np.log2(x/M)[:120], "o-", ms=2, lw=0.5, color="0.25")
ax[0].set_xlabel("$t$"); ax[0].set_ylabel(r"$\log_2(x_t/M)$")
ax[0].set_title("(a) Christoffel cycle $r=15601$, first 120 elements", fontsize=8)
ax[0].set_ylim(-0.05, 1.05)
# panel (b)
rr = np.array([row[0] for row in rows]); br = np.array([row[6] for row in rows])
ax[1].semilogx(rr, br, "o-", ms=2.5, color="0.1", lw=0.8, label=r"$B_r/r$, Theorem 20.2")
for y, lab, ls in [((1+float(TH))/3, "run bound $(1+\\theta)/3$", "--"),
                   (0.75, "Hercher, Theorem 27", ":"),
                   (1/(2*np.log(2)), "$1/(2\\log 2)$", "-.")]:
    ax[1].axhline(y, color="0.45", ls=ls, lw=0.9, label=lab)
ax[1].set_ylim(0.70, 1.01); ax[1].set_xlabel("rung $r$"); ax[1].set_ylabel("constant")
ax[1].set_title("(b) Constant in $M < c\\,F(r)$", fontsize=9)
ax[1].legend(fontsize=6.5, frameon=False, loc="upper right")
fig.tight_layout()
fig.savefig("../fig_ch20_cycle.pdf"); fig.savefig("../fig_ch20_cycle.png", dpi=150)
print("saved")
