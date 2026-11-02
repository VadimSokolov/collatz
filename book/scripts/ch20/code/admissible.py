"""Uniform sampling of admissible walks (K_t <= floor(t theta + delta_c), k_t >= 1,
K_0 = 0, K_r = s) by exact dynamic programming over (t, K), and the statistic
Phi(w) = sum_t 2^{W_t}, W_t = K_t - t theta, which equals 3 x_0 (2^{delta_c} - 1)
for the rational cycle of w. Prints quantiles of Phi/B_r(delta_c)."""
import numpy as np, mpmath as mp, sys, math
mp.mp.dps = 40
TH = mp.log(3, 2)
def run(r, s, nsamp=2000, seed=1):
    d = s - r*TH
    cap = np.array([int(mp.floor(t*TH + d)) for t in range(r+1)])   # K_t <= cap[t]
    Kmax = s
    # logN[t, K]: log number of admissible prefixes ending at K at time t
    logN = np.full((r+1, Kmax+1), -np.inf, dtype=np.float32)
    logN[0, 0] = 0.0
    for t in range(1, r+1):
        prev = logN[t-1]
        m = prev.max()
        if not np.isfinite(m): break
        c = np.cumsum(np.exp(prev - m))            # sum over K' <= K
        with np.errstate(divide="ignore"):
            cur = np.log(np.concatenate(([0.0], c[:-1]))) + m   # K' < K
        cur[cap[t]+1:] = -np.inf
        logN[t] = cur
    logN[r, :s] = -np.inf; logN[r, s+1:] = -np.inf
    assert np.isfinite(logN[r, s])
    # count of admissible words
    rng = np.random.default_rng(seed)
    thf = float(TH)
    Bt = sum(min(float(2**d), 2.0**(float(cap[t]) - t*thf)) for t in range(r))
    vals = []
    for _ in range(nsamp):
        K = s; phi = 0.0
        for t in range(r, 0, -1):
            w = logN[t-1, :K]
            p = np.exp(w - w.max()); p /= p.sum()
            Kp = rng.choice(K, p=p)
            phi += 2.0**(Kp - (t-1)*thf)
            K = Kp
        vals.append(phi/Bt)
    vals = np.array(vals)
    q = np.quantile(vals, [0.01, 0.5, 0.99])
    print(f"r={r} s={s} log10(#admissible)={logN[r,s]/math.log(10):.2f}  Phi/B_r: "
          f"mean={vals.mean():.4f} q01={q[0]:.4f} median={q[1]:.4f} q99={q[2]:.4f} max={vals.max():.4f}")
    return vals
if __name__ == "__main__" and len(sys.argv) == 1:
    for r, s in [(41, 65), (306, 485), (2301, 3647)]:
        run(r, s)

def figure():
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "serif", "font.size": 9})
    fig, ax = plt.subplots(figsize=(6.4, 2.5))
    for (r, s), c in zip([(41, 65), (306, 485), (2301, 3647)], ["0.15", "0.45", "0.7"]):
        v = run(r, s)
        np.save(f"../data/phi_ratio_r{r}.npy", v)
        ax.hist(v, bins=60, histtype="step", color=c, density=True, label=f"$r={r}$")
    ax.axvline(0.75, color="0.3", ls="--", lw=0.9)
    ax.text(0.76, ax.get_ylim()[1]*0.8, "needed at $r^*$, $X=2^{71}$", fontsize=7.5)
    ax.set_xlim(0, 1); ax.set_xlabel(r"$\Phi(w)/B_r$ for a uniform random admissible word $w$")
    ax.set_ylabel("density"); ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig("../fig_ch20_phi.pdf"); fig.savefig("../fig_ch20_phi.png", dpi=150)
if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "fig":
    figure()
