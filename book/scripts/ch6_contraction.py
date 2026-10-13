import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
# Per step contraction factor versus a = alpha*lambda for GD, NAG(0.585), NAG(0.9), SMM (Chapter 6 figure)
import sys; sys.path.insert(0, BOOK + '/scripts')
from bookstyle import *
import numpy as np
from scipy.optimize import brentq
beta = np.log(1.5) / np.log(2)
def rN(a, b): return max(abs(np.roots([1, -(1 - a) * (1 + b), (1 - a) * b])))
def rS(a, b=beta): return np.sqrt(abs((1 - a) * (1 - (1 + b) * a)))
A = np.logspace(-3, 0, 900)
gd = 1 - A; n5 = np.array([rN(a, beta) for a in A]); n9 = np.array([rN(a, 0.9) for a in A]); s = np.array([rS(a) for a in A])
for a in [0.01, 0.04, 0.1, 0.25, 0.5, 1.0]:
    print(a, round(1 - a, 4), round(rN(a, beta), 4), round(rN(a, 0.9), 4), round(rS(a), 4))
c5 = brentq(lambda a: rS(a) - rN(a, beta), 0.15, 0.4); print("crossover SMM=NAG585 at a=", c5)
c9 = brentq(lambda a: rS(a) - rN(a, 0.9), 0.05, 0.2); print("crossover SMM=NAG9 at a=", c9)
cg = 2 / (2 + beta); print("crossover SMM=GD at a=", cg)
fig, ax = plt.subplots(figsize=(WIDTH, 2.9))
best = np.minimum.reduce([gd, n5, n9, s])
ax.fill_between(A, 0, 1.02, where=(s <= best + 1e-12) & (A > 0.2), color=PAL['gold'], alpha=0.13, lw=0)
ax.plot(A, gd, color=PAL['green'], lw=1.8, label='gradient descent')
ax.plot(A, n5, color=PAL['violet'], lw=1.6, ls='--', label=r'Nesterov, $\beta=0.585$')
ax.plot(A, n9, color=PAL['sky'], lw=1.6, ls='-.', label=r'Nesterov, $\beta=0.9$')
ax.plot(A, s, color=PAL['coral'], lw=2.2, label=r'SMM, $\beta=\log_2 3-1$')
for x, lab in [(1 / (1 + beta), r'$a=1/(1+\beta)$')]:
    ax.annotate(lab, (x, 0.0), xytext=(0.16, 0.42), fontsize=8, color=PAL['ink'], ha='center',
                arrowprops=dict(arrowstyle='->', color=PAL['grey'], lw=0.7))
ax.text(np.sqrt(c5 * cg), 0.93, 'SMM fastest', ha='center', fontsize=7.5, color=PAL['ink'])
ax.set_xscale('log'); ax.set_xlabel(r'normalized step $a=\alpha\lambda$')
ax.set_ylabel('per step contraction factor'); ax.set_ylim(0, 1.02); ax.set_xlim(1e-3, 1)
ax.legend(loc='lower left')
fig.savefig(BOOK + '/figures/ch6_contraction.pdf')
