"""Figures and numbers for Chapter 9: Syracuse law mod 9, fine scale TV,
running kbar along long single orbits."""
import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys as _s; _s.path.insert(0, BOOK + '/scripts'); from bookstyle import *
OUT = BOOK + '/figures/'
LT = np.log2(3)
def syr(x):
    m = 3*x+1; k = (m & -m).bit_length()-1
    return m >> k, k

def exact_law(a):
    M = 3**a; units = [r for r in range(M) if r % 3]; order = 2*M//3; inv2 = pow(2, -1, M)
    id_ = {r: i for i, r in enumerate(units)}
    P = np.zeros((len(units),)*2)
    for r in units:
        y = (3*r+1) % M
        for k in range(1, order+1):
            y = y*inv2 % M; P[id_[y], id_[r]] += 2.0**-k/(1-2.0**-order)
    v = np.zeros(len(units)); v[0] = 1
    for _ in range(a): v = P @ v
    return dict(zip(units, v))
pis = {a: exact_law(a) for a in range(1, 8)}
res9 = [1, 2, 4, 5, 7, 8]
law9 = np.array([pis[2][r] for r in res9])

# empirical: large random starts (first 100 strides after 2 burn in strides)
rng = np.random.default_rng(3)
cnt_big = np.zeros(9)
for i in range(3000):
    x = int.from_bytes(rng.bytes(64), 'big') | (1 << 511) | 1
    for j in range(102):
        x, k = syr(x)
        if j >= 2: cnt_big[x % 9] += 1
f_big = cnt_big[res9]/cnt_big.sum()
# empirical: small starting values, all odd iterates after the start, orbits run to 1
cnt_small = np.zeros(9); seen = set(); n_orb = n_end5 = 0
for x0 in range(3, 20002, 2):
    x = last = x0; n_orb += 1
    while x != 1:
        last = x
        x, k = syr(x)
        if x != 1: cnt_small[x % 9] += 1; seen.add(x)
    n_end5 += (last == 5)
f_small = cnt_small[res9]/cnt_small.sum()
print('law9*63', np.round(law9*63, 3)); print('big', np.round(f_big, 4)); print('small', np.round(f_small, 4))
print('TV big vs law', 0.5*np.abs(f_big-law9).sum(), ' TV small vs law', 0.5*np.abs(f_small-law9).sum())
print('mass on 0,3,6 (big, small):', cnt_big[[0,3,6]].sum(), cnt_small[[0,3,6]].sum())
print('small starts: pooled odd iterates %d, distinct values %d, share of orbits whose last stride starts at 5 %.4f'
      % (cnt_small.sum(), len(seen), n_end5/n_orb))

tv_lift = []; tv_unif = []
for a in range(2, 8):
    Mp = 3**(a-1)
    tv_lift.append(0.5*sum(abs(p - pis[a-1][u % Mp]/3) for u, p in pis[a].items()))
    tv_unif.append(0.5*sum(abs(p - 1/len(pis[a])) for p in pis[a].values()))

fig, axes = plt.subplots(1, 2, figsize=(6.4, 2.6))
ax = axes[0]; w = 0.27; xs = np.arange(6); plt.rcParams['hatch.linewidth'] = 0.8
ax.bar(xs - w, law9, w, color='0.3', label='Syracuse law (exact)')
ax.bar(xs, f_big, w, facecolor=PAL['teal'], edgecolor='white', linewidth=0, hatch='///', label='large random starts')  # hatched: teal and dark grey print as one gray
ax.bar(xs + w, f_small, w, color=PAL['gold'], label='orbits of $3\\leq x_0\\leq 20001$')
ax.axhline(1/6, color='k', ls=':', lw=0.8, label='uniform on units')
ax.set_xticks(xs); ax.set_xticklabels([str(r) for r in res9])
ax.set_xlabel('residue of odd iterate mod 9'); ax.set_ylabel('frequency')
ax.legend(frameon=False, loc='upper left')
ax.text(0.62, 0.98, '(a)', transform=ax.transAxes, ha='right', va='top'); ax.set_ylim(0, 0.47)
ax = axes[1]; aa = np.arange(2, 8)
ax.plot(aa, tv_lift, 'o-', ms=3.5, lw=0.9, color=PAL['teal'], label='to lift of law mod $3^{a-1}$')
ax.plot(aa, tv_unif, 's--', ms=3.5, lw=0.9, color=PAL['coral'], label='to uniform on units mod $3^a$')
ax.set_xlabel('$a$'); ax.set_ylabel('total variation'); ax.set_ylim(0, 0.45)
ax.legend(frameon=False, loc='center right')
ax.text(0.98, 0.98, '(b)', transform=ax.transAxes, ha='right', va='top')
plt.tight_layout(); plt.savefig(OUT + 'ch9_syrlaw.pdf', bbox_inches='tight'); plt.close()

# running kbar along long single orbits
fig, axes = plt.subplots(1, 2, figsize=(6.4, 2.6))
rng = np.random.default_rng(11)
cols = [PAL['teal'], PAL['green'], PAL['violet']]
lss = ['-', '-', (0, (4, 1.5, 1, 1.5))]  # teal and violet print as one gray: the 2^4000 orbit is dash-dot
allk = None
for i, bits in enumerate([1000, 2000, 4000]):
    x = int.from_bytes(rng.bytes(bits//8), 'big') | (1 << (bits-1)) | 1
    ks = []
    while x != 1:
        x, k = syr(x); ks.append(k)
    ks = np.array(ks); run = np.cumsum(ks)/np.arange(1, len(ks)+1)
    print('orbit from 2^%d: strides %d, kbar final %.4f, min running kbar after 100 strides %.4f, max %.4f'
          % (bits, len(ks), run[-1], run[100:].min(), run[100:].max()))
    print('   predicted strides', bits*np.log(2)/abs(np.log(0.75)))
    axes[0].plot(np.arange(1, len(ks)+1), run, color=cols[i], lw=0.8, ls=lss[i], label='start near $2^{%d}$' % bits)
    if bits == 4000: allk = ks
axes[0].axhline(2, color='k', lw=0.6, ls='--'); axes[0].axhline(LT, color=PAL['coral'], lw=0.8, ls=':')
axes[0].text(1.1, LT+0.02, r'$\log_2 3$', color=PAL['coral'], ha='left', va='bottom', fontsize=8)  # left end, clear of the legend
axes[0].set_ylim(1.3, 2.6); axes[0].set_xscale('log')
axes[0].set_xlabel('stride $N$'); axes[0].set_ylabel(r'$\bar k_N$')
axes[0].legend(frameon=False, loc='lower right', handlelength=2.5)
for l in axes[0].get_legend().get_lines(): l.set_linewidth(1.6)
kk = np.arange(1, 11)
obs = np.array([(allk == k).mean() for k in kk])
print('k freq along 2^4000 orbit:', np.round(obs, 4), ' n strides', len(allk), ' mean', allk.mean(), ' var', allk.var())
axes[1].bar(kk - 0.18, obs, 0.36, color=PAL['violet'], label='one orbit, start near $2^{4000}$')
axes[1].bar(kk + 0.18, 0.5**kk, 0.36, color='0.6', label='Geometric$(1/2)$')
axes[1].set_yscale('log'); axes[1].set_xticks(kk)
axes[1].set_xlabel('$k_j = v_2(3y_{j-1}+1)$'); axes[1].set_ylabel('frequency')
axes[1].legend(frameon=False, loc='upper right')
plt.tight_layout(); plt.savefig(OUT + 'ch9_kbar.pdf', bbox_inches='tight'); plt.close()
