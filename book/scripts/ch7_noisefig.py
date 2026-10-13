import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
# Chapter 7: the noise robustness experiment (kappa = 25, d = 25, T = 300, R = 30; script ch7_experiments.py, whose
# every step rows use the Nesterov form) drawn as a chart, with the stationary predictions of the exact theory
# (see ch7_variance.py for the predicted ratios).
import sys; sys.path.insert(0, BOOK + '/scripts')
from bookstyle import *
import numpy as np, json
ex = json.load(open(BOOK + '/scripts/ch7_experiments.json'))['noise']
sig = np.array(ex['sigma'][1:])
key = {'GD': 'GD', 'SMM': 'SMM', 'every, 0.5': 'NAG 0.5', 'every, 0.585': 'NAG 0.585', 'every, 0.667': 'NAG 0.667',
       'every, 0.9': 'NAG 0.9'}
obs = {k: [m for m, s in ex['final'][v][1:]] for k, v in key.items()}
pred = {'GD': 1.0, 'every, 0.5': 1.4858, 'every, 0.585': 1.6521, 'SMM': 1.2124, 'every, 0.667': 1.8671, 'every, 0.9': 3.3602}
g0 = 0.30550
sty = {'GD': (PAL['green'], 'o', 'gradient descent'), 'SMM': (PAL['coral'], 'D', r'SMM$_{\mathrm{det}}$, $\beta=0.585$'),
       'every, 0.5': (PAL['sky'], 's', r'every step, $\beta=0.5$'), 'every, 0.585': (PAL['violet'], 'v', r'every step, $\beta=0.585$'),
       'every, 0.667': (PAL['plum'], 'P', r'every step, $\beta=0.667$'), 'every, 0.9': (PAL['gold'], '^', r'every step, $\beta=0.9$')}
order = ['GD', 'SMM', 'every, 0.5', 'every, 0.585', 'every, 0.667', 'every, 0.9']
# Grayscale: violet, plum and sky print alike, so each series has its own marker shape (sizes balance the visual weight).
msz = {'o': 4.2, 'D': 4.0, 's': 4.0, 'v': 4.8, 'P': 5.0, '^': 4.8}
fig, axs = plt.subplots(1, 2, figsize=(WIDTH, 2.55), gridspec_kw=dict(wspace=0.42))
ax = axs[0]; ss = np.geomspace(0.04, 0.5, 50)
for k in order:
    c, m, lab = sty[k]
    ax.plot(ss, g0 * pred[k] * ss ** 2, color=c, lw=1.0, alpha=0.8)
    ax.plot(sig, obs[k], m, color=c, ms=msz[m], mec='white', mew=0.5, label=lab)
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel(r'noise level $\sigma$'); ax.set_ylabel('final objective gap')
ax.set_xticks([0.05, 0.1, 0.2, 0.4]); ax.set_xticklabels(['0.05', '0.1', '0.2', '0.4']); ax.minorticks_off()
ax.text(-0.3, 1.02, '(a)', transform=ax.transAxes, va='bottom')
ax = axs[1]
for j, k in enumerate(order):
    c, m, lab = sty[k]; r = np.array(obs[k]) / np.array(obs['GD'])
    x = j + np.linspace(-0.24, 0.24, 4)
    ax.plot([j - 0.36, j + 0.36], [pred[k]] * 2, color=c, lw=2.4, solid_capstyle='butt')
    ax.plot(x, r, m, color=c, ms=0.86 * msz[m], mec='white', mew=0.4)
ax.set_xticks(range(6)); ax.set_xticklabels(['GD', 'SMM', '0.5', '0.585', '0.667', '0.9'], fontsize=7, rotation=40)
ax.set_xlabel(r'method (every step: $\beta$)'); ax.set_ylabel('gap relative to GD'); ax.set_ylim(0.8, 3.8)
ax.grid(axis='x', visible=False)
ax.text(-0.3, 1.02, '(b)', transform=ax.transAxes, va='bottom')
h, l = axs[0].get_legend_handles_labels()
fig.legend(h, l, loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=3, fontsize=7, columnspacing=1.0, handletextpad=0.3)
fig.savefig(BOOK + '/figures/ch7_noise.pdf')
