import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
os.makedirs(BOOK + '/scripts/data', exist_ok=True)  # intermediate results that later scripts load
# Single orbits in practice (chap:divergence, problem O2): residues of Syracuse iterates modulo 3^a
# along long orbits started near 2^N, compared with the stationary law pi_a of the model
# chain y -> (3y+1) 2^{-K} mod 3^a, K ~ Geometric(1/2), and with i.i.d. sampling from pi_a.
import sys; sys.path.insert(0, BOOK + '/scripts')
from bookstyle import *
import numpy as np, random, json
def pi_a(a):
    M = 3 ** a; cls = [r for r in range(M) if r % 3]
    idx = {r: i for i, r in enumerate(cls)}; n = len(cls)
    P = np.zeros((n, n)); inv2 = pow(2, -1, M)
    for r in cls:
        v = (3 * r + 1) % M
        for K in range(1, 2 * M + 1):            # weights 2^-K, period of 2 mod 3^a divides 2*3^(a-1)
            v = (v * inv2) % M
            P[idx[v], idx[r]] += 2.0 ** (-K)
    P /= P.sum(axis=0, keepdims=True)
    w, V = np.linalg.eig(P); s = np.real(V[:, np.argmin(abs(w - 1))]); s /= s.sum()
    return cls, idx, s
def orbit_residues(y, a):
    M = 3 ** a; out = []
    while y > 1:
        y = 3 * y + 1
        while y % 2 == 0: y //= 2
        out.append(y % M)
    return out
random.seed(7); A = 5
cls, idx, pi = pi_a(A); n = len(cls)
res = {}
for N in (1000, 3000, 10000):
    curves = []
    for rep in range(8):
        y = (1 << N) + random.getrandbits(N - 1) | 1
        r = orbit_residues(y, A)
        cnt = np.zeros(n); tv = []; checkpoints = np.unique(np.geomspace(10, len(r), 40).astype(int))
        j0 = 0
        for c in checkpoints:
            for v in r[j0:c]: cnt[idx[v]] += 1
            j0 = c; tv.append(0.5 * np.abs(cnt / c - pi).sum())
        curves.append((checkpoints.tolist(), tv, len(r)))
    res[N] = curves
    print(N, 'strides', [c[2] for c in curves], 'final TV', [round(c[1][-1], 4) for c in curves], flush=True)
# iid reference
rng = np.random.default_rng(1); ref = []
for c in np.unique(np.geomspace(10, 40000, 40).astype(int)):
    tvs = [0.5 * np.abs(np.bincount(rng.choice(n, c, p=pi), minlength=n) / c - pi).sum() for _ in range(200)]
    ref.append((c, np.mean(tvs)))
ref = np.array(ref)
json.dump({'res': {str(k): v for k, v in res.items()}, 'ref': ref.tolist()}, open(BOOK + '/scripts/data/ch9_orbits.json', 'w'))
fig, axs = plt.subplots(1, 2, figsize=(WIDTH, 2.3), gridspec_kw=dict(wspace=0.35))
ax = axs[0]
cols = {1000: PAL['sky'], 3000: PAL['violet'], 10000: PAL['coral']}
lss = {1000: '-', 3000: (0, (4, 1.5, 1, 1.5)), 10000: (0, (1, 1.6))}  # the three hues print as two grays: separate the start heights by dashes
for N, curves in res.items():
    for i, (cp, tv, L) in enumerate(curves):
        ax.loglog(cp, tv, color=cols[N], lw=0.8, ls=lss[N], alpha=0.8, label=f'start near $2^{{{N}}}$' if i == 0 else None)
ax.loglog(ref[:, 0], ref[:, 1], color=PAL['ink'], lw=1.6, ls='--', label='i.i.d. from $\\pi_5$')
ax.set_xlabel('strides $j$'); ax.set_ylabel('total variation to $\\pi_5$'); ax.legend(fontsize=6)
ax.set_title('(a)', loc='left', fontsize=8)
ax = axs[1]
y = (1 << 10000) + random.getrandbits(9999) | 1
r = np.array(orbit_residues(y, 2)); cls2, idx2, pi2 = pi_a(2)
W = 400; nb = len(r) // W
H = np.array([[np.mean(r[b * W:(b + 1) * W] == c) for b in range(nb)] for c in cls2])
im = ax.imshow(H / pi2[:, None], aspect='auto', cmap=CMAP_SEQ, vmin=0.6, vmax=1.4,
               extent=[0, nb * W, len(cls2) - 0.5, -0.5])
ax.set_yticks(range(len(cls2))); ax.set_yticklabels([str(c) for c in cls2])
ax.set_xlabel('stride $j$'); ax.set_ylabel('residue mod 9'); ax.grid(False)
cb = fig.colorbar(im, ax=ax, fraction=0.05, ticks=[0.6, 0.8, 1.0, 1.2, 1.4]); cb.set_label('frequency / $\\pi_2$', fontsize=7)
cb.ax.axhline(1.0, color='w', lw=1.2)  # sequential map (lighter is more frequent); the line marks frequency = pi_2
ax.set_title('(b)', loc='left', fontsize=8)
fig.savefig(BOOK + '/figures/ch9_orbits.pdf'); plt.close(fig)
print('pi_2', dict(zip(cls2, np.round(pi2 * 63, 3))), 'iid TV at 1e4', np.interp(1e4, ref[:, 0], ref[:, 1]))
