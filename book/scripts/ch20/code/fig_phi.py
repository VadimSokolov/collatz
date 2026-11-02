"""Figure 20.2 from the saved samples of admissible.py."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "serif", "font.size": 9})
fig, ax = plt.subplots(figsize=(6.4, 2.5))
for r, c in zip([41, 306, 2301], ["0.1", "0.45", "0.7"]):
    v = np.load(f"../data/phi_ratio_r{r}.npy")
    bins = np.logspace(-3, 0, 70)
    ax.hist(v, bins=bins, histtype="step", color=c, label=f"$r={r}$")
ax.axvline(0.7513, color="0.3", ls="--", lw=0.9)
ax.text(0.73, 60, "needed at $r^*$\nfor $M\\geq 2^{71}$", fontsize=7.5, ha="right")
ax.set_xscale("log"); ax.set_xlim(1e-3, 1)
ax.set_xlabel(r"$\Phi(w)/B_r$ for a uniform random admissible word $w$")
ax.set_ylabel("count of 2000"); ax.legend(frameon=False, fontsize=8, loc="upper left")
fig.tight_layout(); fig.savefig("../fig_ch20_phi.pdf"); fig.savefig("../fig_ch20_phi.png", dpi=150)
