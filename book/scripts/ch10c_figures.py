import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
import sys; sys.path.insert(0, BOOK + '/scripts')
from bookstyle import *
import numpy as np, json, math
FIG = BOOK + '/figures/'
D = math.log(2) + (math.log(2)/math.log(3))*math.log(math.log(2)/math.log(3)) + (1 - math.log(2)/math.log(3))*math.log(1 - math.log(2)/math.log(3))  # log 2 - H(1/log2 3)
pf = json.load(open(BOOK + '/scripts/data/ch10c_perfect.json')); Z = np.array(pf['Z']); lv = np.array(pf['levels']); ta = np.array(pf['taus'])
mc = json.load(open(BOOK + '/scripts/data/ch10c_mc.json')); ms = np.array(mc['ms']); ex = np.array(mc['res']['exact']); exactZ = np.array(mc['exactZ'])
fig, axs = plt.subplots(1, 3, figsize=(WIDTH, 2.3), gridspec_kw=dict(wspace=0.75))
ax = axs[0]; k = np.arange(1, 401)
ax.semilogy(k, Z[1:], color=PAL['teal'], lw=1.5, label=r'$Z(m)=P(\tau>m)$ exact')
ax.semilogy(k, np.exp(-D * k), color=PAL['grey'], lw=0.8, ls=':', label=r'$e^{-Dm}$')
A = np.vstack([np.ones(301), -np.log(k[99:])]).T; c = np.linalg.lstsq(A, np.log(Z[100:401]) + D * k[99:], rcond=None)[0]
ax.semilogy(k[20:], np.exp(c[0]) * k[20:] ** (-c[1]) * np.exp(-D * k[20:]), color=PAL['coral'], lw=0.9, ls='--', label=r'$C\,m^{-%.3f}e^{-Dm}$' % c[1])
ax.set_xlabel('depth $m$'); ax.set_ylabel('volume'); ax.legend(fontsize=5.5, loc='lower left'); ax.set_title('(a)', loc='left', fontsize=8)
ax = axs[1]
ax.plot(np.arange(3000), lv[:3000], color=PAL['violet'], lw=0.5)
ax.set_xlabel('sweep'); ax.set_ylabel('level $k$'); ax.set_title('(b) ideal weight chain', loc='left', fontsize=8)
ax = axs[2]
mm = np.arange(0, 601)
ax.semilogy(mm, exactZ, color=PAL['ink'], lw=1.2, label='exact')
cols = {'naive': PAL['grey'], 'nested': PAL['teal'], 'split': PAL['coral'], 'vertical': PAL['gold']}
mks = {'naive': ('s', 3.3), 'nested': ('o', 3.5), 'split': ('^', 4.2), 'vertical': ('D', 3.1)}  # one shape per estimator: the light colours print as similar grays
for name in ('naive', 'nested', 'split', 'vertical'):
    est = np.array(mc['res'][name]); med = np.median(est, axis=0); med[med <= 0] = np.nan
    ax.semilogy(ms, med, mks[name][0], color=cols[name], ms=mks[name][1], mec=PAL['ink'], mew=0.4, label=name, alpha=0.9)
ax.set_xlabel(r'$m$'); ax.set_ylabel(r'$P(\sigma_\infty>m)$'); ax.set_ylim(1e-9, 1.5); ax.legend(fontsize=5.5, loc='lower left')
ax.set_title('(c) beyond the horizon', loc='left', fontsize=8)
fig.savefig(FIG + 'ch10c_montecarlo.pdf'); plt.close(fig)
print('fit', c)
