import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
# Chapter 7: figures of the experiments (condition number sweep, frequency ablation, logistic regression, minibatch
# gradients), drawn from ch7_experiments.json, which ch7_experiments.py writes.
import sys; sys.path.insert(0, BOOK + '/scripts')
from bookstyle import *
import numpy as np, json

r = json.load(open(BOOK + '/scripts/ch7_experiments.json'))
FIG = BOOK + '/figures/'
BC = {'0.5': PAL['sky'], '0.585': PAL['violet'], '0.667': PAL['teal'], '0.9': PAL['gold']}
SMM = r'SMM$_{\mathrm{det}}$'

# ---------------------------------------------------------------- condition number sweep
k = r['ksweep']; K = np.array(k['kappa']); C = k['counts']; cap = k['T'] + 1
fig, ax = plt.subplots(figsize=(WIDTH, 2.9))
ax.axhline(cap, color=PAL['grey'], lw=0.8, ls=':')
ax.text(2.05, cap * 1.08, r'cap $T+1$', fontsize=7, color=PAL['ink'])
# Grayscale: the colours of beta = 0.585 and 0.667 print alike, so the marker fill alternates (open, filled, open, filled) and the size grows with beta.
KF = {'0.5': (False, 3.6), '0.585': (True, 3.8), '0.667': (False, 4.4), '0.9': (True, 4.6)}
for nm in ['0.5', '0.585', '0.667', '0.9']:
    fl, sz = KF[nm]; mf = BC[nm] if fl else 'white'
    ax.plot(K, [c[0] for c in C['HB ' + nm]], color=BC[nm], marker='s', ms=sz, mfc=mf, mew=0.9, lw=1.1,
            label=rf'heavy ball, $\beta={nm}$')
    ax.plot(K, [c[0] for c in C['NAG ' + nm]], color=BC[nm], marker='^', ms=sz + 0.5, mfc=mf, mew=0.9, lw=1.1, ls='--',
            label=rf'Nesterov, $\beta={nm}$')
ax.plot(K, [c[0] for c in C['GD']], color=PAL['ink'], marker='o', ms=4.0, lw=1.8, label='gradient descent')
ax.plot(K, [c[0] for c in C['SMM']], color=PAL['coral'], marker='D', ms=4.2, lw=2.2, label=SMM + r', $\beta=0.585$')
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(1.8, 115); ax.set_ylim(3.5, 500)
ax.set_xticks(K); ax.set_xticklabels([str(x) for x in K]); ax.minorticks_off()
ax.set_xlabel(r'condition number $\kappa$'); ax.set_ylabel(r'iterations to a gap below $10^{-4}$')
h, l = ax.get_legend_handles_labels(); o = [8, 9, 0, 1, 2, 3, 4, 5, 6, 7]
ax.legend([h[i] for i in o], [l[i] for i in o], loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=3,
          fontsize=7, columnspacing=1.0, handlelength=2.2)
fig.savefig(FIG + 'ch7_ksweep.pdf')

# ---------------------------------------------------------------- frequency ablation
fig, ax = plt.subplots(figsize=(WIDTH, 2.6))
ax.axhline(cap, color=PAL['grey'], lw=0.8, ls=':')
ax.text(2.05, cap * 1.08, r'cap $T+1$', fontsize=7, color=PAL['ink'])
for key, lab, c, m, ls in [('NAG 0.585', r'every step, $\beta=0.585$', PAL['violet'], 'o', '-'),
                           ('SMM', 'odd steps, ' + SMM, PAL['coral'], 'D', '-'),
                           ('random', r'random steps, probability $\frac{1}{2}$', PAL['gold'], 's', '--'),
                           ('GD', 'never, gradient descent', PAL['ink'], 'x', ':')]:
    v = np.array(C[key])
    ax.errorbar(K, v[:, 0], yerr=v[:, 1], color=c, marker=m, ms=4.2, lw=1.6, ls=ls, capsize=2, label=lab)
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(1.8, 115); ax.set_ylim(3.5, 500)
ax.set_xticks(K); ax.set_xticklabels([str(x) for x in K]); ax.minorticks_off()
ax.set_xlabel(r'condition number $\kappa$'); ax.set_ylabel(r'iterations to a gap below $10^{-4}$')
ax.legend(loc='lower right', fontsize=7)
fig.savefig(FIG + 'ch7_ablation.pdf')

# ---------------------------------------------------------------- logistic regression
lg = r['logistic']; cv = lg['curves']
fig, ax = plt.subplots(figsize=(WIDTH, 2.7))
t = np.arange(len(cv['GD']))
for key, lab, c, ls, lw in [('GD', 'gradient descent', PAL['ink'], '-', 1.6),
                            ('HB tuned', 'heavy ball, tuned', PAL['sky'], '--', 1.4),
                            ('NAG tuned', 'Nesterov, tuned', PAL['gold'], '-.', 1.4),
                            ('SMM', SMM, PAL['coral'], '-', 2.0),
                            ('SMM state', r'SMM$_{\mathrm{state}}$', PAL['violet'], '-', 1.4)]:
    ax.plot(t, cv[key], color=c, ls=ls, lw=lw, label=lab)
ax.axhline(lg['tol'], color=PAL['grey'], lw=0.8, ls=':')
ax.set_yscale('log'); ax.set_xlim(0, len(t) - 1); ax.set_ylim(1e-17, 2)
ax.set_xlabel(r'iteration $t$'); ax.set_ylabel(r'objective gap $f(w_t)-f^\ast$')
ax.legend(loc='upper right', fontsize=7)
fig.savefig(FIG + 'ch7_logistic.pdf')

# ---------------------------------------------------------------- minibatch gradients
sg = r['sgd']
fig, axs = plt.subplots(1, 3, figsize=(WIDTH, 2.25), sharey=True, gridspec_kw=dict(wspace=0.12))
sty = [('GD', 'gradient descent', PAL['ink'], '-', 1.2), ('SMM', SMM, PAL['coral'], '-', 1.4),
       ('NAG 0.585', r'Nesterov, $\beta=0.585$', PAL['violet'], '-', 1.0),
       ('NAG tuned', 'Nesterov, tuned', PAL['sky'], '--', 1.0), ('NAG 0.9', r'Nesterov, $\beta=0.9$', PAL['gold'], '-', 1.0)]
for j, (ax, B) in enumerate(zip(axs, sg['batch'])):
    for key, lab, c, ls, lw in sty:
        y = np.array(sg['curves'][key][j]); ax.plot(np.arange(len(y)), y, color=c, ls=ls, lw=lw, label=lab)
    ax.axhline(sg['tol'], color=PAL['grey'], lw=0.8, ls=':')
    ax.set_yscale('log'); ax.set_xlim(0, sg['T']); ax.set_ylim(3e-5, 1)
    ax.set_xticks([0, 200, 400, 600]); ax.tick_params(labelsize=7)
    ax.text(0.95, 0.95, rf'$|B|={B}$', transform=ax.transAxes, ha='right', va='top', fontsize=7.5)
    ax.set_xlabel(r'iteration $t$', fontsize=8)
axs[0].set_ylabel('objective gap', fontsize=8)
h, l = axs[0].get_legend_handles_labels()
fig.legend(h, l, loc='upper center', bbox_to_anchor=(0.5, -0.02), ncol=3, fontsize=7, columnspacing=1.0)
fig.savefig(FIG + 'ch7_sgd.pdf')
print('figures written')
