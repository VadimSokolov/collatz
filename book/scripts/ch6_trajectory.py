import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
# Chapter 6 trajectory figure: tilted two dimensional quadratic, kappa = 8.
# f(x) = x^T H x / 2, H has eigenvalues L = 8 (eigenvector at 55 degrees) and mu = 1; x0 = (1.3, 0.9); alpha = 1/L.
# Methods: GD, SMM_det (Algorithm smmdet), SMM_state exactly as Algorithm smmstate (f_ref = 0, rho_thr = 0.356),
# and every step Nesterov with beta = log2(3) - 1.
import sys; sys.path.insert(0, BOOK + '/scripts')
from bookstyle import *
import numpy as np
b = np.log(1.5) / np.log(2); L = 8.0; T = 28
th = np.deg2rad(55.0); v = np.array([np.cos(th), np.sin(th)]); w = np.array([-v[1], v[0]])
H = L * np.outer(v, v) + np.outer(w, w)
f = lambda x: 0.5 * x @ H @ x; g = lambda x: H @ x
x0 = np.array([1.3, 0.9])

def run(kind):
    x = x0.copy(); xp = x0.copy(); gp = max(f(x0), 1e-12); X = [x.copy()]; inj = []
    for t in range(T):
        if kind == 'gd': do = False
        elif kind == 'det': do = (t % 2 == 1)
        elif kind == 'nag': do = True
        else:
            gt = max(f(x), 1e-12); do = (t == 0 or gt / gp > 0.356); gp = gt
        y = x + b * (x - xp) if do else x
        if do: inj.append(t)
        xp, x = x, y - g(y) / L; X.append(x.copy())
    return np.array(X), inj

R = {k: run(k) for k in ['gd', 'det', 'state', 'nag']}
for k, (X, inj) in R.items():
    print(k, 'f(x_T) = %.3e' % f(X[-1]), 'injections at', inj)
fig, axs = plt.subplots(2, 1, figsize=(WIDTH, 4.6), gridspec_kw=dict(height_ratios=[1.0, 0.95], hspace=0.38))
ax = axs[0]
X0, X1, Y0, Y1 = -0.06, 0.46, -0.31, 0.07
gx, gy = np.meshgrid(np.linspace(X0, X1, 400), np.linspace(Y0, Y1, 400))
Z = 0.5 * (H[0, 0] * gx ** 2 + 2 * H[0, 1] * gx * gy + H[1, 1] * gy ** 2)
ax.contourf(gx, gy, np.log(Z + 1e-3), levels=18, cmap='Blues_r', alpha=0.35)
ax.contour(gx, gy, Z, levels=np.geomspace(2e-5, 0.4, 16), colors=PAL['teal'], linewidths=0.35, alpha=0.6)
col = dict(gd=PAL['green'], det=PAL['coral'], state=PAL['violet'], nag=PAL['sky'])
lab = dict(gd='GD', det=r'SMM$_{\mathrm{det}}$', state=r'SMM$_{\mathrm{state}}$', nag=r'NAG, $\beta=0.585$')
for k in ['gd', 'state', 'det']:
    X, inj = R[k]
    ax.plot(X[:, 0], X[:, 1], color=col[k], lw=1.3 if k == 'gd' else 1.0, ls='-.' if k == 'gd' else '-', alpha=0.9)  # GD dash-dot: grayscale-safe
    if k == 'det':
        idx = np.arange(1, T + 1)
        after_inj = np.array([i - 1 in inj for i in idx])
        ax.scatter(X[idx[~after_inj], 0], X[idx[~after_inj], 1], s=16, facecolor='white', edgecolor=col[k], lw=0.8, zorder=4)
        ax.scatter(X[idx[after_inj], 0], X[idx[after_inj], 1], s=18, marker='D', color=col[k], lw=0, zorder=4)
    else:
        ax.scatter(X[1:, 0], X[1:, 1], s=6, color=col[k], zorder=3)
ax.annotate(r'first step from $x_0=(1.3,\,0.9)$', xy=(0.447, -0.245), xytext=(0.0, -0.265), fontsize=7, color=PAL['ink'], va='center', arrowprops=dict(arrowstyle='->', color=PAL['ink'], lw=0.6))
ax.plot(0, 0, '*', color=PAL['gold'], ms=9, mec=PAL['ink'], mew=0.5, zorder=5)
ax.set_xlim(X0, X1); ax.set_ylim(Y0, Y1); ax.set_aspect('equal'); ax.grid(False)
ax.set_xlabel(r'$x_1$'); ax.set_ylabel(r'$x_2$')
ax = axs[1]
for k in ['gd', 'det', 'state', 'nag']:
    X, _ = R[k]; ax.semilogy(range(T + 1), [f(x) for x in X], color=col[k], lw=1.6, label=lab[k],
                             ls={'nag': '--', 'gd': '-.'}.get(k, '-'))
ax.set_xlabel(r'iteration $t$'); ax.set_ylabel(r'$f(x_t)$')
ax.legend(loc='lower left', fontsize=7.5)
h = [plt.Line2D([], [], color=col[k], lw=1.4, ls='-.' if k == 'gd' else '-', label=lab[k]) for k in ['gd', 'state', 'det']]
h += [plt.Line2D([], [], marker='o', ls='', mfc='white', mec=PAL['coral'], ms=4, label='after descent step'), plt.Line2D([], [], marker='D', ls='', color=PAL['coral'], ms=4, label='after injection step')]
axs[0].legend(handles=h, loc='upper left', bbox_to_anchor=(1.03, 1.0), fontsize=7, ncol=1)
fig.savefig(BOOK + '/figures/ch6_trajectory.pdf')
