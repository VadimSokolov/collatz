"""Figures for the chapter on progress on the open problems.
ch16_o4.pdf     : (a) third cumulant against E log n: census, sampled blocks, model and fit; (b) skewness.
ch16_smm.pdf    : (a) worst case of Switching MM over the main formula, two branches; (b) noise floors.
ch16_erdos.pdf  : (a) Erdos 125 density just below 3^m against frac(m log_4 3); (b) single orbit statistics."""
import sys, os, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bookstyle import *
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); FIG = os.path.join(HERE, '..', 'figures'); DATA = os.path.join(HERE, 'data')
L2 = math.log(2)

# ---------------------------------------------------------------- O4
o4 = json.load(open(os.path.join(DATA, 'ch16_o4.json'))); pr = json.load(open(os.path.join(DATA, 'ch16_o4_predict.json')))
rows = np.load(os.path.join(DATA, 'ch3_decomp.npz'))['rows']
ks = rows[:, 0]; eu = rows[:, 2]; k3 = rows[:, 7]; k2 = rows[:, 9]
f2, f3 = o4['fit']['k2'], o4['fit']['k3']
fig, axs = plt.subplots(1, 2, figsize=(WIDTH, 2.3), gridspec_kw=dict(wspace=0.42))
ax = axs[0]; sel = ks >= 10
uu = np.linspace(6, 36, 100)
ax.plot(uu, f3['model_slope'] * uu / 1e3, color=PAL['grey'], lw=0.9, ls='--', label='model')
ax.plot(uu, (f3['slope'] * uu - f3['D']) / 1e3, color=PAL['teal'], lw=0.9, label='model minus constant')
ax.plot(eu[sel], k3[sel] / 1e3, 'o', color=PAL['ink'], ms=2.5, label='census, $n<2^{30}$')
pu = [(r['k'] + 0.5) * L2 for r in pr]
ax.plot(pu, [r['k3'] / 1e3 for r in pr], 's', color=PAL['coral'], ms=3, label='sampled, $k\\leq48$')
ax.set_xlabel(r'$\mathbb{E}\log n$'); ax.set_ylabel(r'third cumulant ($10^3$)'); ax.legend(fontsize=5.5, loc='upper left')
ax.set_title('(a)', loc='left', fontsize=8)
ax = axs[1]; kk = np.arange(12, 70)
u2 = (kk + 0.5) * L2
fitsk = (f3['slope'] * u2 - f3['D']) / (f2['slope'] * u2 - f2['D']) ** 1.5
modsk = (f3['model_slope'] * u2) / (f2['model_slope'] * u2) ** 1.5
ax.plot(kk, modsk, color=PAL['grey'], lw=0.9, ls='--', label='model')
ax.plot(kk, fitsk, color=PAL['teal'], lw=0.9, label='fit from $k\\leq29$')
ax.plot(ks[ks >= 12], (k3 / k2 ** 1.5)[ks >= 12], 'o', color=PAL['ink'], ms=2.5, label='census')
ax.errorbar([r['k'] for r in pr], [r['skew'] for r in pr], yerr=[2 * r['skew_se'] for r in pr], fmt='s', color=PAL['coral'], ms=3, lw=0.8, label='sampled')
ax.set_ylim(0.4, 1.5); ax.set_xlabel('block $k$, $n\\in[2^k,2^{k+1})$'); ax.set_ylabel('skewness of stopping time')
ax.legend(fontsize=5.5, loc='upper right'); ax.set_title('(b)', loc='left', fontsize=8)
fig.savefig(os.path.join(FIG, 'ch16_o4.pdf')); plt.close(fig)

# ---------------------------------------------------------------- SMM: two branches and noise floors
fc = json.load(open(os.path.join(DATA, 'ch16_fullconj.json'))); s2 = json.load(open(os.path.join(DATA, 'ch16_s2.json')))
fig, axs = plt.subplots(1, 2, figsize=(WIDTH, 2.3), gridspec_kw=dict(wspace=0.45))
ax = axs[0]; bb = np.linspace(0, 0.995, 300)
cols = {2: PAL['coral'], 4: PAL['violet'], 6: PAL['teal'], 8: PAL['gold'], 10: PAL['sky']}
for N, col in cols.items():
    main = 1 / ((4 + 2 * bb) * N + 2); over = bb ** 2 / (4 * N - 2 + 2 * (N - 2) * bb)
    ax.plot(bb, np.maximum(main, over) / main, color=col, lw=1.6, label=f'$N={N}$')
    pts = [r for r in fc if r[0] == N]
    ax.plot([r[1] for r in pts], [r[2] * ((4 + 2 * r[1]) * N + 2) for r in pts], 'o', color=col, ms=2.6)
ax.axvline(math.log(1.5) / math.log(2), color=PAL['ink'], lw=0.6); ax.text(0.57, 1.5, 'Collatz', fontsize=6, rotation=90, ha='right', va='center')
ax.set_xlabel(r'$\beta$'); ax.set_ylabel('worst case / main branch'); ax.set_ylim(0.97, 1.9)
ax.legend(fontsize=5.5, loc='upper left'); ax.set_title('(a)', loc='left', fontsize=8)
ax = axs[1]; meth = ['gd', 'smm', 'nag585', 'nag9', 'hb9']
lab = ['GD', 'SMM', 'NAG\n0.585', 'NAG\n0.9', 'HB\n0.9']
x = np.arange(len(meth)); w = 0.2
for j, (b, col) in enumerate((('16', PAL['coral']), ('64', PAL['teal']))):
    R = s2['res'][b]; gd_p, gd_s = R['pred']['gd'], R['sim']['gd'][0]
    ax.bar(x + (2 * j - 1.5) * w, [R['pred'][m] / gd_p for m in meth], w, color=col, alpha=0.35, label=f'linear theory, $b={b}$')
    ax.bar(x + (2 * j - 0.5) * w, [R['sim'][m][0] / gd_s for m in meth], w, color=col, alpha=0.95, label=f'simulation, $b={b}$')
ax.axhline(5 / 3, color=PAL['ink'], lw=0.7, ls=':'); ax.text(1.35, 5 / 3 * 1.08, '$5/3$', fontsize=6)
ax.set_yscale('log'); ax.set_xticks(x); ax.set_xticklabels(lab, fontsize=6); ax.set_ylabel('noise floor / gradient descent')
ax.legend(fontsize=5, loc='upper left'); ax.set_title('(b)', loc='left', fontsize=8)
fig.savefig(os.path.join(FIG, 'ch16_smm.pdf')); plt.close(fig)

# ---------------------------------------------------------------- Erdos 125 and single orbits
er = json.load(open(os.path.join(DATA, 'ch16_erdos125_17.json'))); p3 = json.load(open(os.path.join(DATA, 'ch16_p3.json')))
fig, axs = plt.subplots(1, 2, figsize=(WIDTH, 2.3), gridspec_kw=dict(wspace=0.42))
ax = axs[0]; L43 = math.log(3) / math.log(4)
for r in er['pow3']:
    f = (r['m'] * L43) % 1.0
    ax.plot(f, r['below'], 'o', color=PAL['teal'] if r['m'] <= 16 else PAL['coral'], ms=3.2)
    ax.text(f + 0.012, r['below'] + 0.002, str(r['m']), fontsize=5.5)
ax.axhline(1015 / 1458, color=PAL['grey'], lw=0.7, ls=':'); ax.text(0.02, 1015 / 1458 + 0.004, 'Hasler and Melfi bound', fontsize=5.5)
ax.plot([], [], 'o', color=PAL['teal'], ms=3, label=r'$m\leq16$ (to $4^{13}$)'); ax.plot([], [], 'o', color=PAL['coral'], ms=3, label=r'$17\leq m\leq21$ (to $4^{17}$)')
ax.set_xlim(0, 1); ax.set_ylim(0.68, 0.9); ax.set_xlabel(r'fractional part of $m\log_43$'); ax.set_ylabel(r'density at $3^m-1$')
ax.legend(fontsize=5.5, loc='lower center', bbox_to_anchor=(0.5, 0.08)); ax.set_title('(a)', loc='left', fontsize=8)
ax = axs[1]
for i, rec in enumerate(p3):
    col = PAL['coral'] if rec['bits'] == 300000 else PAL['teal']
    a = [5, 6, 7, 8]; ratio = [rec['tv'][str(q)]['tv'] / rec['tv'][str(q)]['iid_mean'] for q in a]
    ax.plot(np.array(a) + 0.06 * (i - 2), ratio, 'o', color=col, ms=3)
ax.axhline(1, color=PAL['ink'], lw=0.7)
ax.plot([], [], 'o', color=PAL['teal'], ms=3, label=r'four orbits from $2^{100000}$'); ax.plot([], [], 'o', color=PAL['coral'], ms=3, label=r'one orbit from $2^{300000}$')
ax.set_xticks([5, 6, 7, 8]); ax.set_xlabel(r'modulus $3^a$'); ax.set_ylabel('distance to $\\pi_a$ / i.i.d. value')
ax.set_ylim(0.8, 1.2); ax.legend(fontsize=5.5, loc='upper right'); ax.set_title('(b)', loc='left', fontsize=8)
fig.savefig(os.path.join(FIG, 'ch16_erdos.pdf')); plt.close(fig)
print('ok')
