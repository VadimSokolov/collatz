import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
import sys; sys.path.insert(0, BOOK + '/scripts')
from bookstyle import *
import numpy as np, json, math
FIG = BOOK + '/figures/'
a = json.load(open(BOOK + '/scripts/data/ch10c_acceptance.json'))
cl = json.load(open(BOOK + '/scripts/data/ch10c_clusters.json'))
ms = np.array(a['ms']); Z = np.array(a['Z'])
fig, axs = plt.subplots(1, 2, figsize=(WIDTH, 2.2), gridspec_kw=dict(wspace=0.45, width_ratios=[1.25, 1]))
ax = axs[0]
ax.semilogy(ms, Z, color=PAL['ink'], lw=1.2, label=r'$Z_\sigma(m)$, rejection sampling')
ax.semilogy(ms, a['acc'], 'o-', color=PAL['teal'], ms=3, lw=1, label='all bits')
lo = np.array(a['acc_lo']); hi = np.array(a['acc_hi']); hi[hi <= 0] = np.nan
ax.semilogy(ms, lo, 's-', color=PAL['coral'], ms=2.5, lw=0.9, label='bits $0$ to $12$')
ax.semilogy(ms, hi, '^-', color=PAL['gold'], ms=2.5, lw=0.9, label='bits $13$ to $25$')
ax.set_ylim(1e-7, 2); ax.set_xlabel('level $m$'); ax.set_ylabel('acceptance of a bit flip')
ax.legend(fontsize=5.5, loc='lower left'); ax.set_title('(a)', loc='left', fontsize=8)
ax = axs[1]
cols = {'300': PAL['teal'], '350': PAL['coral'], '400': PAL['violet'], '440': PAL['gold']}
for m in ('300', '350', '400', '440'):
    s = np.array(cl[m]['sizes'], float); r = np.arange(1, len(s) + 1)
    ax.loglog(r, s, 'o-', color=cols[m], ms=2.2, lw=0.8, label=r'$m=%s$: %d classes' % (m, cl[m]['clusters']))
ax.set_xlabel('rank of class'); ax.set_ylabel('size'); ax.legend(fontsize=5.5, loc='upper right')
ax.set_title('(b)', loc='left', fontsize=8)
fig.savefig(FIG + 'ch10c_acceptance.pdf'); plt.close(fig)
