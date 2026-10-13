"""Chapter 3: stopping time slope and CLT figure (dash free labels).
Left: total stopping time of the accelerated map T versus log n, n <= 1e4,
with the model slope 1/|delta*| and the least squares slope through the origin.
Right: standardized residuals (tau - L/|d|)/sqrt(s2 L/|d|^3), pooled over
4000 random n in each of [1e3,2e3), ..., [1e7,2e7), against N(0,1)."""
import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys as _s; _s.path.insert(0, BOOK + '/scripts'); from bookstyle import *
d = abs(0.5*np.log(0.75)); s2 = 0.25*np.log(3)**2
def T(n): return n >> 1 if n % 2 == 0 else (3*n+1) >> 1
def tau(n):
    s = 0
    while n != 1: n = T(n); s += 1
    return s
ns = np.arange(2, 10001); st = np.array([tau(int(n)) for n in ns]); L = np.log(ns)
g = (L @ st)/(L @ L)
rng = np.random.default_rng(1)
Z = []
for lo in [1e3, 1e4, 1e5, 1e6, 1e7]:
    sample = rng.integers(int(lo), int(2*lo), 4000)
    s = np.array([tau(int(n)) for n in sample]); l = np.log(sample)
    Z.append((s - l/d)/np.sqrt(s2*l/d**3))
Z = np.concatenate(Z)
print("slope %.4f; pooled z: mean %.3f var %.3f skew %.3f; frac |z|>2 %.4f (normal 0.0455)" %
      (g, Z.mean(), Z.var(), ((Z-Z.mean())**3).mean()/Z.std()**3, (np.abs(Z) > 2).mean()))
fig, axes = plt.subplots(1, 2, figsize=(6.4, 2.6))
ax = axes[0]
ax.plot(L, st, '.', ms=1.0, color=PAL['teal'], alpha=0.6, rasterized=True)
t = np.linspace(0, L.max(), 50)
ax.plot(t, t/d, 'k--', lw=0.9, zorder=3, label=r'model slope $1/|\delta^*| = 6.952$')  # drawn above the fitted line, which nearly coincides with it
ax.plot(t, g*t, color=PAL['coral'], lw=1.6, zorder=2, label=r'fitted slope $%.3f$' % g)
ax.set_xlabel(r'$\log n$'); ax.set_ylabel('total stopping time of $T$')
ax.legend(frameon=False, loc='upper left')
ax = axes[1]
ax.hist(Z, bins=50, density=True, color=PAL['teal'], alpha=0.6, label='orbits, $n$ in five ranges')
u = np.linspace(-4, 5, 200); ax.plot(u, np.exp(-u*u/2)/np.sqrt(2*np.pi), 'k', lw=0.9, label='standard normal')
ax.set_xlabel('standardized total stopping time'); ax.set_ylabel('density')
ax.legend(frameon=False, loc='upper right')
plt.tight_layout(); plt.savefig(BOOK + '/figures/ch3_stopping.pdf', bbox_inches='tight', dpi=300)
