import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
# Period two cyclic steps (h1, h2) on spectrum [mu, L]: worst case two step factor R, slow mode
# noise floor factor N = 2(h1^2+h2^2)/(h1+h2)^2 at matched speed, and the Pareto frontier.
import numpy as np, sympy as sp
from scipy.optimize import brentq
b0 = np.log(1.5) / np.log(2)

def R_grid(h1, h2, mu, L, n=200001):
    lam = np.linspace(mu, L, n); return np.abs((1 - h1 * lam) * (1 - h2 * lam)).max()
def N_of(h1, h2): return 2 * (h1 ** 2 + h2 ** 2) / (h1 + h2) ** 2
def G_of_r(r): return (1 + r) ** 2 / (4 * r)
def Rstar(k): return ((k - 1) / (k + 1)) ** 2            # best constant step 2/(L+mu), two steps
def Gcheb(k): return 2 * (k + 1) ** 2 / (k ** 2 + 6 * k + 1)
def r_of_G(G): return 2 * G - 1 + 2 * np.sqrt(G * (G - 1))
def Rfront(N, k): return 1 - 4 * k / ((2 - N) * (k + 1) ** 2)

# 1. random check that every pair lies on or above the frontier, with equality at the predicted scale
rng = np.random.default_rng(11); worst = np.inf
for k in [2.0, 5.0, 10.0, 30.0, 100.0]:
    mu, L = 1.0, k
    for _ in range(400):
        h1, h2 = rng.uniform(0.01, 4.0 / L, 2)
        R = R_grid(h1, h2, mu, L, 20001); N = min(N_of(h1, h2), (3 * k * k + 2 * k + 3) / (2 * (k + 1) ** 2))
        worst = min(worst, R - Rfront(N, k))
print("min over random pairs of R - Rfront(N) (should be >= 0 up to grid error):", worst)
for k in [5.0, 100.0]:
    for r in [1.0, 1 + b0, 3.0]:
        mu, L = 1.0, k; h2 = (1 + r) / (L + mu); h1 = h2 / r
        print(f"kappa={k} r={r:.3f}: R grid {R_grid(h1,h2,mu,L):.6f} frontier {Rfront(N_of(h1,h2),k):.6f}  G={G_of_r(r):.4f}")

# 2. symbolic checks
k, r = sp.symbols('kappa r', positive=True)
Gc = 2 * (k + 1) ** 2 / (k ** 2 + 6 * k + 1)
print("1 - Gcheb*(1-R*) == Chebyshev value:", sp.simplify(1 - Gc * 4 * k / (k + 1) ** 2 - (k - 1) ** 2 / (k ** 2 + 6 * k + 1)) == 0)
print("Ncheb = 2 - 1/Gcheb:", sp.factor(2 - 1 / Gc))
print("N(r) = 2 - 1/G(r):", sp.simplify(2 * (1 + r ** 2) / (1 + r) ** 2 - (2 - 4 * r / (1 + r) ** 2)) == 0)

# 3. kappa at which the Chebyshev ratio equals 1 + beta
GS = G_of_r(1 + b0)
k0 = brentq(lambda kk: Gcheb(kk) - GS, 1.0001, 10)
kk = sp.symbols('kk', positive=True); Gs = sp.nsimplify(GS)
print("G_SMM =", GS, " N_SMM =", 2 - 1 / GS, " kappa0 =", k0)
# closed form: Gcheb(kappa) = G  <=>  2(k+1)^2 = G (k^2+6k+1)
Gsym = sp.symbols('G', positive=True)
sol = sp.solve(sp.Eq(2 * (kk + 1) ** 2, Gsym * (kk ** 2 + 6 * kk + 1)), kk)
print("roots:", [sp.simplify(s) for s in sol], [float(s.subs(Gsym, GS)) for s in sol])

# 4. table
print("\nkappa, r_cheb, N_cheb, G_cheb, per step factor Cheb, SMM(1/L): R, G, per step ; SMM optimally scaled alpha*L, per step")
tab = []
for kap in [2, 5, 10, 25, 100, 1000]:
    mu, L = 1.0, float(kap)
    Gc_ = Gcheb(kap); rc = r_of_G(Gc_); Nc = 2 - 1 / Gc_
    Rc = (kap - 1) ** 2 / (kap ** 2 + 6 * kap + 1)
    Rs = R_grid(1 / L, (1 + b0) / L, mu, L); Gs_ = (1 - Rs) / (1 - Rstar(kap))
    Rgd = R_grid(1 / L, 1 / L, mu, L); Ggd = (1 - Rgd) / (1 - Rstar(kap))
    alpha_star = (2 + b0) / ((1 + b0) * (L + mu))
    Ropt = R_grid(alpha_star, (1 + b0) * alpha_star, mu, L)
    onfront = GS <= Gc_
    tab.append((kap, rc, Nc, Gc_, np.sqrt(Rc), Rs, Gs_, np.sqrt(Rs), Ggd, alpha_star * L, np.sqrt(Ropt), onfront))
    print(f"{kap:5d} r_cheb={rc:.3f} N_cheb={Nc:.4f} G_cheb={Gc_:.4f} cheb/step={np.sqrt(Rc):.4f} | SMM 1/L: G={Gs_:.4f} step={np.sqrt(Rs):.4f} (GD 1/L G={Ggd:.4f}) | SMM opt: alpha*L={alpha_star*L:.4f} step={np.sqrt(Ropt):.4f} frontier? {onfront}")
print("limits: r_cheb ->", r_of_G(2.0), " N_cheb -> 1.5, alpha* L -> (2+b)/(1+b) =", (2 + b0) / (1 + b0))
print("SMM exponent gain at optimal scale G_SMM =", GS, "; at alpha=1/L, kappa->inf: (2+b)/4 =", (2 + b0) / 4)

# 5. figure
import sys; sys.path.insert(0, BOOK + '/scripts'); from bookstyle import *
fig, axs = plt.subplots(2, 1, figsize=(WIDTH, WIDTH * 1.18))
ax = axs[0]
kaps = [2, 5, 10, 100]; cols = [PAL['sky'], PAL['teal'], PAL['violet'], PAL['plum']]
Nends = [1.0] + [2 - 1 / Gcheb(x) for x in kaps]
NN = np.linspace(1, 1.5, 400); ax.plot(NN, 1 / (2 - NN), color=PAL['grey'], lw=1.0, ls=':', zorder=1)
for j, (kap, c) in enumerate(zip(kaps, cols)):
    n = np.linspace(Nends[j], Nends[j + 1], 100)
    ax.plot(n, 1 / (2 - n), color=c, lw=3.0, solid_capstyle='butt', zorder=2)
    ax.plot([Nends[j + 1]], [Gcheb(kap)], 'o', color=c, ms=6.5, mec='white', mew=0.9, zorder=4,
            label=rf'Chebyshev end, $\kappa={kap}$')
ax.plot([1.5], [2.0], 'o', color=PAL['ink'], ms=6, mfc='white', zorder=4)
ax.annotate(r'Chebyshev, $\kappa\to\infty$', (1.5, 2.0), xytext=(1.33, 1.88), fontsize=7.5, color=PAL['ink'],
            ha='center', va='center', arrowprops=dict(arrowstyle='->', color=PAL['grey'], lw=0.7))
ax.plot([2 - 1 / GS], [GS], marker='*', color=PAL['gold'], ms=14, mec=PAL['ink'], mew=0.6, zorder=5)
ax.annotate('SMM ratio, best scale', (2 - 1 / GS, GS), xytext=(1.13, 1.45), fontsize=7.5, color=PAL['ink'],
            arrowprops=dict(arrowstyle='->', color=PAL['grey'], lw=0.7))
ax.plot([1.0], [1.0], 's', color=PAL['green'], ms=7, mec='white', zorder=5)
ax.annotate(r'GD, step $2/(L+\mu)$', (1.0, 1.0), xytext=(1.06, 0.86), fontsize=7.5, color=PAL['ink'],
            arrowprops=dict(arrowstyle='->', color=PAL['grey'], lw=0.7))
Rs = R_grid(1 / 100, (1 + b0) / 100, 1, 100); g_s = (1 - Rs) / (1 - Rstar(100))
Rg = R_grid(1 / 100, 1 / 100, 1, 100); g_g = (1 - Rg) / (1 - Rstar(100))
ax.plot([1.0], [g_g], 's', color=PAL['green'], ms=6, mfc='white', zorder=5)
ax.plot([2 - 1 / GS], [g_s], marker='*', color=PAL['gold'], ms=12, mfc='white', mec=PAL['gold'], zorder=5)
ax.annotate(r'SMM, $\alpha=1/L$, $\kappa=100$', (2 - 1 / GS, g_s), xytext=(1.12, 0.70), fontsize=7.5,
            color=PAL['ink'], arrowprops=dict(arrowstyle='->', color=PAL['grey'], lw=0.7))
ax.annotate(r'GD, $\alpha=1/L$, $\kappa=100$', (1.0, g_g), xytext=(1.12, 0.50), fontsize=7.5,
            color=PAL['ink'], arrowprops=dict(arrowstyle='->', color=PAL['grey'], lw=0.7))
ax.set_xlim(0.98, 1.53); ax.set_ylim(0.4, 2.15)
ax.set_xlabel(r'noise floor at matched speed $N$')
ax.set_ylabel(r'rate gain $G$')
ax.legend(loc='upper left', fontsize=7)
ax.text(-0.13, 1.0, '(a)', transform=ax.transAxes, va='top')
ax = axs[1]
K = np.logspace(0.0001, 3, 400); rc = [r_of_G(Gcheb(x)) for x in K]
ax.fill_between(K, 1, rc, color=PAL['violet'], alpha=0.12, lw=0)
ax.plot(K, rc, color=PAL['violet'], lw=2.2, label=r'Chebyshev ratio $r_{\mathrm{ch}}(\kappa)$')
ax.axhline((np.sqrt(2) + 1) ** 2, color=PAL['violet'], lw=0.8, ls=':')
ax.text(250, 5.05, r'$(\sqrt{2}+1)^2$', fontsize=8, color=PAL['violet'])
ax.axhline(1 + b0, color=PAL['gold'], lw=2.2, label=r'SMM ratio $1+\beta$')
ax.axvline(k0, color=PAL['grey'], lw=0.8, ls='--')
ax.text(k0 * 1.12, 1.15, rf'$\kappa_0={k0:.2f}$', fontsize=7.5, color=PAL['ink'])
ax.text(40, 3.6, 'frontier ratios', fontsize=7.5, color=PAL['violet'], ha='center')
ax.set_xscale('log'); ax.set_ylim(1, 6.6); ax.set_xlim(1, 1000)
ax.set_xlabel(r'condition number $\kappa$'); ax.set_ylabel(r'step ratio $h_2/h_1$')
ax.legend(loc='lower right', fontsize=7.5, bbox_to_anchor=(1.0, 0.12))
ax.text(-0.13, 1.0, '(b)', transform=ax.transAxes, va='top')
fig.tight_layout(); fig.savefig(BOOK + '/figures/ch7_pareto.pdf')
