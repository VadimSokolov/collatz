"""Divergence chapter (chap:divergence), single orbits in practice: checks behind
the text of sec:ch9:single and rem:ch9:null (no figure).
(1) The 24 orbits of ch9_orbits.py (same seed): strides to reach 1, final total
    variation to pi_5, and each running total variation compared with the spread
    of independent samples of the same size from pi_5.
(2) The model chain on the units mod 3^5 (independent Geometric(1/2) halving
    exponents, started from pi_5): mean running total variation after 2400, 7000
    and 24000 strides, against independent draws from pi_5.
(3) Chi squared statistics of the residues mod 9 along each orbit of the chapter
    (the 24 orbits above and the three orbits of fig:ch9:kbar) against the
    uniform law on the six units and on all nine residues, and against pi_2;
    calibration of the chi squared test against pi_2 under the model chain.
(4) The running total variation of the 24 orbits compared with the spread of the
    model chain at the same checkpoints (its successive residues are dependent,
    so its spread exceeds that of independent draws)."""
import random
import numpy as np
from scipy.stats import chi2


def pi_a(a):
    """Stationary law of the model chain on the units mod 3^a (as in ch9_orbits.py)."""
    M = 3 ** a; cls = [r for r in range(M) if r % 3]
    idx = {r: i for i, r in enumerate(cls)}; n = len(cls)
    P = np.zeros((n, n)); inv2 = pow(2, -1, M)
    for r in cls:
        v = (3 * r + 1) % M
        for K in range(1, 2 * M + 1):
            v = (v * inv2) % M
            P[idx[v], idx[r]] += 2.0 ** (-K)
    P /= P.sum(axis=0, keepdims=True)
    w, V = np.linalg.eig(P); s = np.real(V[:, np.argmin(abs(w - 1))]); s /= s.sum()
    return cls, idx, s


def orbit_odd_iterates(y):
    out = []
    while y > 1:
        y = 3 * y + 1
        while y % 2 == 0: y //= 2
        out.append(y)
    return out


def tv(p, q): return 0.5 * np.abs(p - q).sum()


A = 5; M = 3 ** A
cls, idx, pi = pi_a(A); n = len(cls)
lut = np.full(M, -1); lut[cls] = np.arange(n)
rng = np.random.default_rng(2026)

# (1) the 24 orbits of ch9_orbits.py, same random stream
random.seed(7)
orbits = {}
for N in (1000, 3000, 10000):
    orbits[N] = []
    for rep in range(8):
        y = (1 << N) + random.getrandbits(N - 1) | 1
        orbits[N].append(orbit_odd_iterates(y))
REPS = 400
worst = []
print('(1) orbits of ch9_orbits.py: residues mod 3^5 against pi_5')
for N, group in orbits.items():
    lens = [len(o) for o in group]; finals = []
    for o in group:
        r = lut[np.array([v % M for v in o])]
        cps = np.unique(np.geomspace(10, len(r), 40).astype(int))
        cnt = np.zeros(n); j0 = 0; zs = []; pct = []; over99 = 0
        for c in cps:
            np.add.at(cnt, r[j0:c], 1); j0 = c
            t = tv(cnt / c, pi)
            sims = 0.5 * np.abs(rng.multinomial(c, pi, size=REPS) / c - pi).sum(axis=1)
            if c >= 100:
                zs.append((t - sims.mean()) / sims.std()); pct.append((sims < t).mean())
                over99 += t > np.quantile(sims, 0.99)
        finals.append(t)
        worst.append((max(zs), max(pct), over99, len(zs), N))
    c_mean = int(np.mean(lens))
    sims = 0.5 * np.abs(rng.multinomial(c_mean, pi, size=2000) / c_mean - pi).sum(axis=1)
    print('  start near 2^%d: strides %d to %d, final TV mean %.4f (range %.4f to %.4f);'
          ' i.i.d. at %d draws: mean %.4f, sd %.4f, 99%% quantile %.4f'
          % (N, min(lens), max(lens), np.mean(finals), min(finals), max(finals), c_mean,
             sims.mean(), sims.std(), np.quantile(sims, 0.99)))
w = np.array([x[:4] for x in worst])
print('  over all 24 orbits and checkpoints j >= 100: largest z score %.2f, largest i.i.d.'
      ' percentile %.3f, checkpoints above the i.i.d. 99%% quantile: %d of %d'
      % (w[:, 0].max(), w[:, 1].max(), int(w[:, 2].sum()), int(w[:, 3].sum())))
print('  per orbit (largest z, largest percentile):',
      [(x[4], round(x[0], 2), round(x[1], 3)) for x in worst])

# (2) the model chain, started from pi_5, against i.i.d. draws
R = 1000; checkpoints = [2400, 7000, 24000]
inv2 = pow(2, -1, M); order = 2 * 3 ** (A - 1)
p2 = np.array([pow(inv2, k, M) for k in range(order)])
units = np.array(cls)
y = rng.choice(units, size=R, p=pi)
cnt = np.zeros((R, n)); rows = np.arange(R); out = {}
for j in range(1, max(checkpoints) + 1):
    K = rng.geometric(0.5, size=R)
    y = ((3 * y + 1) * p2[K % order]) % M
    cnt[rows, lut[y]] += 1
    if j in checkpoints:
        t = 0.5 * np.abs(cnt / j - pi).sum(axis=1); out[j] = (t.mean(), t.std())
print('(2) running TV to pi_5, mean (sd) over %d paths' % R)
for c in checkpoints:
    sims = 0.5 * np.abs(rng.multinomial(c, pi, size=4000) / c - pi).sum(axis=1)
    print('  %5d strides: model chain %.4f (%.4f), i.i.d. %.4f (%.4f)'
          % (c, out[c][0], out[c][1], sims.mean(), sims.std()))

# (3) chi squared tests on residues mod 9
cls2, idx2, pi2 = pi_a(2)
rng11 = np.random.default_rng(11)          # the three orbits of fig:ch9:kbar (ch9_figures.py)
kbar_orbits = []
for bits in [1000, 2000, 4000]:
    x = int.from_bytes(rng11.bytes(bits // 8), 'big') | (1 << (bits - 1)) | 1
    kbar_orbits.append(orbit_odd_iterates(x))
allorb = [o for g in orbits.values() for o in g] + kbar_orbits
res = []
for o in allorb:
    c9 = np.bincount(np.array([v % 9 for v in o]), minlength=9).astype(float); L = c9.sum()
    c6 = c9[cls2]
    s6 = ((c6 - L / 6) ** 2 / (L / 6)).sum(); s9 = ((c9 - L / 9) ** 2 / (L / 9)).sum()
    sp = ((c6 - L * pi2) ** 2 / (L * pi2)).sum()
    res.append((L, s6, chi2.sf(s6, 5), s9, chi2.sf(s9, 8), sp, chi2.sf(sp, 5)))
res = np.array(res)
print('(3) chi squared on residues mod 9, %d orbits of %d to %d strides' % (len(res), res[:, 0].min(), res[:, 0].max()))
print('  against uniform on the six units: smallest statistic %.0f, largest p value %.1e'
      % (res[:, 1].min(), res[:, 2].max()))
print('  against uniform on all nine residues: smallest statistic %.0f, largest p value %.1e'
      % (res[:, 3].min(), res[:, 4].max()))
print('  against pi_2: statistics %.1f to %.1f, p values %.3f to %.3f, rejections at 5%%: %d'
      % (res[:, 5].min(), res[:, 5].max(), res[:, 6].min(), res[:, 6].max(), int((res[:, 6] < 0.05).sum())))
# calibration of the chi squared test against pi_2 when the residues follow the model chain mod 9
M2 = 9; inv2 = pow(2, -1, M2); p2 = np.array([pow(inv2, k, M2) for k in range(6)])
lut2 = np.full(M2, -1); lut2[cls2] = np.arange(6)
R = 20000; L = 2000
y = rng.choice(np.array(cls2), size=R, p=pi2); cnt = np.zeros((R, 6)); rows = np.arange(R)
for j in range(L):
    K = rng.geometric(0.5, size=R)
    y = ((3 * y + 1) * p2[K % 6]) % M2
    cnt[rows, lut2[y]] += 1
st = ((cnt - L * pi2) ** 2 / (L * pi2)).sum(axis=1)
iid = rng.multinomial(L, pi2, size=R); sti = ((iid - L * pi2) ** 2 / (L * pi2)).sum(axis=1)
print('  model chain mod 9, %d strides, %d paths: mean statistic against pi_2 %.2f (i.i.d. %.2f, nominal 5),'
      ' rejection rate at nominal 5%% %.3f (i.i.d. %.3f)'
      % (L, R, st.mean(), sti.mean(), (chi2.sf(st, 5) < 0.05).mean(), (chi2.sf(sti, 5) < 0.05).mean()))

# (4) the 24 orbits against the spread of the model chain at the same checkpoints
R = 2000
inv2 = pow(2, -1, M); p2 = np.array([pow(inv2, k, M) for k in range(order)])
orbit_res = [lut[np.array([v % M for v in o])] for g in orbits.values() for o in g]
orbit_cps = [np.unique(np.geomspace(10, len(r), 40).astype(int)) for r in orbit_res]
grid = np.unique(np.concatenate(orbit_cps)); gpos = {c: i for i, c in enumerate(grid)}
TVm = np.zeros((len(grid), R))
y = rng.choice(units, size=R, p=pi); cnt = np.zeros((R, n)); rows = np.arange(R); g = 0
for j in range(1, grid[-1] + 1):
    K = rng.geometric(0.5, size=R)
    y = ((3 * y + 1) * p2[K % order]) % M
    cnt[rows, lut[y]] += 1
    if j == grid[g]:
        TVm[g] = 0.5 * np.abs(cnt / j - pi).sum(axis=1); g += 1
zs = []; pcts = []; over = 0; tot = 0; finals = []
for r, cps in zip(orbit_res, orbit_cps):
    cnt1 = np.zeros(n); j0 = 0
    for c in cps:
        np.add.at(cnt1, r[j0:c], 1); j0 = c
        if c < 100: continue
        t = tv(cnt1 / c, pi); ref = TVm[gpos[c]]
        zs.append((t - ref.mean()) / ref.std()); pcts.append((ref < t).mean())
        over += t > np.quantile(ref, 0.99); tot += 1
    finals.append((ref < t).mean())
print('(4) against the model chain (%d paths) at the same checkpoints j >= 100: largest z score %.2f,'
      ' largest percentile %.4f, checkpoints above its 99%% quantile: %d of %d' % (R, max(zs), max(pcts), over, tot))
print('  percentiles of the final total variation of each orbit in the model chain law:', np.round(sorted(finals), 3).tolist())
