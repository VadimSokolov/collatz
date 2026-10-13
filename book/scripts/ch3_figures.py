"""Figures for Chapter 3: model paths versus orbits, parity ACF, Ville bound."""
import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys as _s; _s.path.insert(0, BOOK + '/scripts'); from bookstyle import *
OUT = BOOK + '/figures/'
LOG2 = np.log(2); ds = 0.5*np.log(0.75)
def T(n): return n >> 1 if n % 2 == 0 else (3*n+1) >> 1
def T0(n): return n >> 1 if n % 2 == 0 else 3*n+1

# ---- Figure: model paths versus real orbits
rng = np.random.default_rng(7)
fig, axes = plt.subplots(1, 2, figsize=(6.4, 2.6), sharey=True)
x0 = 10**12
cols = plt.cm.viridis(np.linspace(0.1, 0.9, 10))
lens_real = []
for i in range(10):
    n = x0 + 2*int(rng.integers(0, 10**6)) + 1
    path = [np.log(n)]
    while n != 1:
        n = T(n); path.append(np.log(n))
    lens_real.append(len(path)-1)
    axes[0].plot(path, color=cols[i], lw=0.8)
lens_mod = []
for i in range(10):
    s = [np.log(x0)]
    while s[-1] > 0:
        s.append(s[-1] + (np.log(1.5) if rng.random() < 0.5 else -LOG2))
    lens_mod.append(len(s)-1)
    axes[1].plot(s, color=cols[i], lw=0.8)
for ax in axes:
    t = np.arange(0, 400)
    ax.plot(t, np.log(x0) + ds*t, 'k--', lw=0.9)
    ax.axhline(0, color='grey', lw=0.5)
    ax.set_xlabel('step $n$ (accelerated map)')
    ax.set_xlim(0, 330); ax.set_ylim(-1.5, 36)
axes[0].set_ylabel(r'$\log x_n$')
axes[0].text(0.97, 0.95, '(a) orbits of $T$', transform=axes[0].transAxes, ha='right', va='top')
axes[1].text(0.97, 0.95, '(b) Bernoulli model', transform=axes[1].transAxes, ha='right', va='top')
plt.tight_layout(); plt.savefig(OUT + 'ch3_paths.pdf', bbox_inches='tight'); plt.close()
print('real lengths', lens_real, np.mean(lens_real)); print('model lengths', lens_mod, np.mean(lens_mod))
print('predicted mean', np.log(x0)/abs(ds), ' sd', np.sqrt(0.301737/abs(ds)**3*np.log(x0)))

# ---- Figure: ACF
rho_large = np.load(BOOK + '/scripts/ch3_acf.npy')
rho_small = [-0.504, 0.2482, -0.153, 0.0822, -0.0032, -0.0139, 0.0019, -0.0031]
H = np.arange(1, 11)
fig, ax = plt.subplots(figsize=(4.6, 2.6))
ax.axhline(0, color='grey', lw=0.5)
ax.plot(H, (-0.5)**H, 'k-', lw=0.9, label=r'renewal model $(-1/2)^h$')
ax.plot(H, rho_large, 'o', ms=4.5, mfc='none', color=PAL['teal'], label='large random starts, steps 50 to 299')
ax.plot(H[:8], rho_small, 's', ms=3.5, color=PAL['coral'], label='small starting values, run to 1')
ax.set_xlabel('lag $h$'); ax.set_ylabel(r'$\rho(h)$'); ax.set_xticks(H)
ax.legend(frameon=False, loc='lower right')
plt.tight_layout(); plt.savefig(OUT + 'ch3_acf.pdf', bbox_inches='tight'); plt.close()

# ---- Figure: Ville bound
ratio = np.load(BOOK + '/scripts/ch3_ratio.npy')
Ks = 2.0**np.arange(1, 13)
emp = np.array([(ratio >= K).mean() for K in Ks])
# model: iid walk, steps log(3/2), -log 2, prob 1/2, sup over all time
W = 1_000_000
S = np.zeros(W); mx = np.zeros(W); alive = np.ones(W, bool)
for t in range(1500):
    step = np.where(rng.random(W) < 0.5, np.log(1.5), -LOG2)
    S = S + step*alive
    mx = np.maximum(mx, S)
    alive &= (S > -40) & (mx < np.log(Ks[-1]) + 1)
mod = np.array([(mx >= np.log(K)).mean() for K in Ks])
print('alive at end', alive.sum())
print('K*emp', np.round(Ks*emp, 3)); print('K*model', np.round(Ks*mod, 3))
fig, ax = plt.subplots(figsize=(4.6, 2.6))
ax.axhline(1, color='k', lw=0.9, label='Ville bound $K\\cdot(1/K)=1$')
ax.plot(Ks, Ks*emp, 'o-', ms=3.5, lw=0.8, color=PAL['teal'], label='orbits of $T$, $2\\leq n\\leq 10^6$')
ax.plot(Ks, Ks*mod, 's--', ms=3.5, lw=0.8, color=PAL['coral'], label='Bernoulli model, $10^6$ walks')
ax.set_xscale('log', base=2); ax.set_ylim(0, 1.15)
ax.set_xlabel('$K$'); ax.set_ylabel(r'$K\cdot \Pr(\sup_n x_n\geq K x_0)$')
ax.legend(frameon=False, loc='lower left')
plt.tight_layout(); plt.savefig(OUT + 'ch3_ville.pdf', bbox_inches='tight'); plt.close()
