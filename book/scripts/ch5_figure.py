"""Figure for Chapter 5: standard Nesterov coefficient versus condition number."""
import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys as _s; _s.path.insert(0, BOOK + '/scripts'); from bookstyle import *
beta = np.log(1.5)/np.log(2)
k = np.geomspace(1, 2000, 600)
c = (np.sqrt(k)-1)/(np.sqrt(k)+1)
fig, ax = plt.subplots(figsize=(4.8, 3.0))
ax.semilogx(k, c, color=PAL['teal'], lw=1.4, label="$(\\sqrt{\\kappa}-1)/(\\sqrt{\\kappa}+1)$")
for b, col, name in [(beta, "#c0392b", "$\\beta=\\log_2 3-1\\approx 0.585$"), (0.9, "#239b56", "$\\beta=0.9$")]:
    kc = ((1+b)/(1-b))**2
    ax.axhline(b, color=col, lw=1.1, ls="--", label=name)
    ax.plot([kc], [b], "o", color=col, ms=4)
    ax.annotate(f"$\\kappa\\approx{kc:.1f}$" if kc < 100 else f"$\\kappa={kc:.0f}$",
                (kc, b), xytext=(6, -12), textcoords="offset points", fontsize=8, color=col)
for kk, lab in [(1.7757, "1.78"), (2.4094, "2.41"), (4.1851, "4.19")]:
    ax.axvline(kk, color="gray", lw=0.5, ls=":")
    ax.text(kk*1.04, 0.02, lab, fontsize=6.5, color="gray", rotation=90, va="bottom")
ax.set_xlabel("condition number $\\kappa$"); ax.set_ylabel("momentum coefficient")
ax.set_ylim(0, 1); ax.set_xlim(1, 2000)
ax.legend(frameon=False, fontsize=7, loc="lower right")
plt.tight_layout()
plt.savefig(BOOK + "/figures/ch5_nesterov_beta.pdf", bbox_inches="tight")
