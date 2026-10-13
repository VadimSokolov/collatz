"""Chapter 8: the window lemma and the refined problem, explored numerically.

Part A  window [L, U] = [(1-(2/3)^m)/(e^Lambda - 1), m/(3 Lambda)] on the trivial
        cycle and on the hypothetical (r, s) of every upper ladder rung up to
        r_*(2^68).  Figure ch8_window.pdf.
Part B  all parity words (or a random sample of 10^6 when there are more than
        3*10^6) with m <= 20 ones and ceil(m theta) <= k <= 2m letters:
        residues c_k(w) mod D, D = 2^k - 3^m, integer hits, closest approach to
        an integer, and a uniform residue model for comparison.
Part C  rotation classes (rational cycles) for m <= 11: position of the class
        minimum in the window, the two means of Proposition means, and the
        extremal words of c_k.
Part D  balanced (Christoffel) words with k = ceil(m theta): least odd, greatest
        odd and greatest element relative to F(m), for the primitive words
        (gcd(m, k) = 1) with 100 <= m <= 2000.  For gcd(m, k) > 1 the word is a
        shorter balanced word repeated and defines the same rational cycle.
Figure ch8_words.pdf.  The right panel divides the least nonzero distance of
c_k(w)/D to an integer, over all words examined, by its expected value when
the residues of those words are independent and uniform on the nonzero
classes modulo D.
"""
import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
import sys; sys.path.insert(0, BOOK + '/scripts'); from bookstyle import *
import numpy as np, math, itertools, mpmath as mp
from math import comb, log
from fractions import Fraction as Fr
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm

mp.mp.dps = 60
TH = mp.log(3)/mp.log(2); L2 = mp.log(2); L3 = mp.log(3)
OUT = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s); OUT.append(s)

def cvalue(word):
    """c_k(w) = sum_j 3^(m-j) 2^(i_j), via c <- 3c + 2^t at each one (eq. ch2:crec)."""
    c = 0
    for t, b in enumerate(word):
        if b: c = 3*c + 2**t
    return c

def window(m, k):
    Lam = k*L2 - m*L3
    lo = (1 - (mp.mpf(2)/3)**m)/(mp.e**Lam - 1)
    return Lam, lo, m/(3*Lam)

# ---------------- Part A ----------------
say('Part A: window checks')
Lam, lo, hi = window(1, 2)
say(f'  trivial cycle m=1 k=2: window [{mp.nstr(lo,6)}, {mp.nstr(hi,6)}], x_min = 1 = c/(2^k-3^m) = {Fr(cvalue([1,0]),1)}')
a = []; y = TH
for _ in range(40):
    ai = int(mp.floor(y)); a.append(ai); y = 1/(y-ai)
q = [1, 0]
for ai in a:
    q.append(ai*q[-1]+q[-2])
Q = lambda k: q[k+2]
rungs = sorted(Q(k-1)+j*Q(k) for k in range(0, 30, 2) for j in range(1, a[k+1]+1))
rungs = [u for u in rungs if u <= 72057431991]
X68 = mp.mpf(2)**68
rows = []
for u in rungs:
    s = int(mp.ceil(u*TH)); Lam, lo, hi = window(u, s)
    rows.append((u, s, lo, hi))
say('  rungs:', len(rows))
say('  r, s, lower L, upper U, (U/L)/(m/3), log2 L, log2 U')
for (u, s, lo, hi) in rows:
    say(f'   {u:>14d} {s:>15d}  {mp.nstr(lo,5):>11s} {mp.nstr(hi,5):>11s}  {mp.nstr((hi/lo)/(mp.mpf(u)/3),8):>11s}'
        f'  {float(mp.log(lo,2)):7.3f} {float(mp.log(hi,2)):7.3f}')
say('  rungs whose lower end exceeds 2^68:', [u for (u, s, lo, hi) in rows if lo > X68])

fig, ax = plt.subplots(figsize=(WIDTH, 3.0))
us = np.array([float(r[0]) for r in rows]); los = np.array([float(r[2]) for r in rows])
his = np.array([float(r[3]) for r in rows])
cmap = plt.get_cmap(CMAP_SEQ)
ax.fill_between(us, los, his, color=PAL['sky'], alpha=0.16, lw=0)
for i, (uu, l, h) in enumerate(zip(us, los, his)):
    ax.plot([uu, uu], [l, h], color=cmap(0.12 + 0.8*i/len(us)), lw=2.0, solid_capstyle='round', alpha=0.95)
ax.plot(us, los, 'v', ms=3.4, color=PAL['violet'], mec=PAL['ink'], mew=0.3,
        label=r'lower end $(1-(2/3)^m)/(e^{\Lambda}-1)$')
ax.plot(us, his, '^', ms=3.4, color=PAL['coral'], mec=PAL['ink'], mew=0.3,
        label=r'upper end $m/(3\Lambda)=F(r)$')
ax.axhline(float(X68), color=PAL['gold'], ls='--', lw=1.2)
ax.text(1.4, float(X68)*4, r'verified range $2^{68}$', fontsize=7, color=PAL['ink'])
ax.annotate(r'$r_*=72\,057\,431\,991$', xy=(us[-1], his[-1]), xytext=(2e5, 3e21), fontsize=7,
            arrowprops=dict(arrowstyle='->', lw=0.6, color=PAL['ink']))
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_ylim(0.3, 1e23)
ax.set_xlabel(r'ladder rung $r=m$, with $k=s=\lceil r\theta\rceil$')
ax.set_ylabel(r'window for $x_{\min}$')
ax.legend(loc='lower right', fontsize=7)
fig.savefig(BOOK + '/figures/ch8_window.pdf'); plt.close(fig)

# ---------------- Part B ----------------
say('\nPart B: parity words, m <= 20, ceil(m theta) <= k <= 2m')
LIMIT = 3_000_000; NSAMP = 1_000_000
rng = np.random.default_rng(20261010)

def enumerate_words(m, k, D):
    """All words of length k with m ones: residues c mod D and c as a float."""
    res = {0: np.zeros(1, dtype=np.int64)}; flt = {0: np.zeros(1)}
    for t in range(k):
        nres, nflt = {}, {}
        p2 = pow(2, t, D); p2f = float(2**t)
        for o in range(0, m+1):
            pr, pf = [], []
            if o in res and o + (k - t - 1) >= m:
                pr.append(res[o]); pf.append(flt[o])
            if (o-1) in res:
                pr.append((3*res[o-1] + p2) % D); pf.append(3*flt[o-1] + p2f)
            if pr:
                nres[o] = np.concatenate(pr); nflt[o] = np.concatenate(pf)
        res, flt = nres, nflt
    return res[m], flt[m]

def sample_words(m, k, D, n):
    pos = np.sort(np.argsort(rng.random((n, k)), axis=1)[:, :m], axis=1)
    pw3 = np.array([pow(3, m-1-j, D) for j in range(m)], dtype=np.int64)
    p2 = np.array([pow(2, t, D) for t in range(k)], dtype=np.int64)
    r = np.zeros(n, dtype=np.int64); f = np.zeros(n)
    for j in range(m):
        r = (r + (pw3[j]*p2[pos[:, j]]) % D) % D
        f += 3.0**(m-1-j)*2.0**pos[:, j]
    return r, f

def expected_min_dist(n, D):
    """Expected least distance to 0 modulo D, divided by D, of n independent
    residues uniform on the nonzero classes modulo the odd number D.  With
    h = (D-1)/2 the distance of one residue is uniform on 1..h, so
    P(min >= j) = ((h-j+1)/h)^n and E[min] = sum_{i=1}^{h} (i/h)^n."""
    h = (D - 1)//2
    if h <= 2*10**7:
        i = np.arange(1, h + 1, dtype=float)
        E = float(np.sum((i/h)**n))
    else:
        E = h/(n + 1) + 0.5          # Euler-Maclaurin; here n is far below h
    return E/D

summary = []; Emodel = 0.0; Emodel_nt = 0.0; Efull = []
for m in range(1, 21):
    k0 = int(mp.ceil(m*TH))
    for k in range(k0, 2*m+1):
        D = 2**k - 3**m
        Lam = k*log(2) - m*log(3); up = m/(3*Lam)
        N = comb(k, m)
        cmin, cmax = 3**m - 2**m, 2**(k-m)*(3**m - 2**m)
        if N <= LIMIT:
            r, f = enumerate_words(m, k, D); mode = 'all'; Neff = N/k
            assert abs(f.min()/cmin - 1) < 1e-12 and abs(f.max()/cmax - 1) < 1e-12
        else:
            r, f = sample_words(m, k, D, NSAMP); mode = 'sample'; Neff = NSAMP
        x = f/D
        nint = int(np.sum(r == 0))
        dist = np.minimum(r, D - r)/D
        nz = dist[r != 0]
        mind = float(nz.min()) if nz.size else float('nan')
        Ez = expected_min_dist(nz.size, D) if nz.size else float('nan')
        rel = mind/Ez
        frac_in = float(np.mean(x < up))
        Emodel += N/(k*D)
        if k < 2*m: Emodel_nt += N/(k*D)
        if k < 2*m and mode == 'all': Efull.append(N/(k*D))
        summary.append((m, k, D, N, mode, nint, mind, 2*Neff*mind, frac_in, up, 2*len(r)*mind, rel, Ez))
say(' m  k               D              words  mode    ints   min dist   2Neff*min  min/E_unif  frac x<U')
for (m, k, D, N, mode, nint, mind, sc, fin, up, scw, rel, Ez) in summary:
    say(f'{m:2d} {k:2d} {D:15d} {N:18d}  {mode:6s} {nint:4d}  {mind:10.3e} {sc:9.3f}  {rel:9.3f}   {fin:.4f}')
hits = [(m, k, nint) for (m, k, D, N, mode, nint, *_r) in summary if nint]
say('integer x hits (m, k, count):', hits)
say('  all hits have k = 2m (words (10)^(m) and (01)^(m), x = 1 and 2):', all(k == 2*m for (m, k, n) in hits))
sc = np.array([s[7] for s in summary if not math.isnan(s[7]) and s[3] >= 50])
say('scaled closeness 2 Neff min dist (N >= 50): n=%d mean %.3f median %.3f min %.4f max %.3f'
    % (len(sc), sc.mean(), np.median(sc), sc.min(), sc.max()))
say('uniform residue model: expected integer classes sum C(k,m)/(k D) over the table = %.3f; excluding k = 2m: %.3f'
    % (Emodel, Emodel_nt))
say('  restricted to fully enumerated (m,k) with k < 2m: %.3f over %d pairs; Poisson P(0) = %.3f'
    % (sum(Efull), len(Efull), math.exp(-sum(Efull))))
scw = np.array([s[10] for s in summary if not math.isnan(s[10]) and s[3] >= 50])
say('per word normalisation 2 N min dist: mean %.3f median %.3f' % (scw.mean(), np.median(scw)))
for lab, sel in [('all', lambda s: True), ('enumerated', lambda s: s[4] == 'all'),
                 ('sampled', lambda s: s[4] == 'sample')]:
    rl = np.array([s[11] for s in summary if not math.isnan(s[11]) and s[3] >= 50 and sel(s)])
    say('min dist / expected under uniform residues (N >= 50), %s: n=%d mean %.3f median %.3f min %.4f max %.3f'
        % (lab, len(rl), rl.mean(), np.median(rl), rl.min(), rl.max()))
rl = [(s[0], s[11]) for s in summary if not math.isnan(s[11]) and s[3] >= 50]
for lo_m, hi_m in [(4, 8), (9, 12), (13, 16), (17, 20), (4, 12), (13, 20)]:
    v = np.array([r for (m, r) in rl if lo_m <= m <= hi_m])
    say('   m in [%d, %d]: mean ratio %.3f, median ratio %.3f over %d shapes'
        % (lo_m, hi_m, v.mean(), np.median(v), len(v)))
ev = [s[12] for s in summary if not math.isnan(s[11]) and s[3] >= 50]
say('   expected closest approach under the model, over these shapes: from %.2e to %.2e' % (max(ev), min(ev)))
fr = np.array([s[8] for s in summary if s[0] >= 5])
say('fraction of words with x(w) < U (m >= 5): min %.4f max %.4f' % (fr.min(), fr.max()))

# ---------------- Part C ----------------
say('\nPart C: rotation classes, m <= 11, all k in [ceil(m theta), 2m]')
cls_stats = []; R11 = {}
for m in range(2, 12):
    k0 = int(mp.ceil(m*TH))
    for k in range(k0, 2*m+1):
        D = 2**k - 3**m; Lam = k*log(2) - m*log(3); up = m/(3*Lam)
        low = (1 - (2/3)**m)/(math.exp(Lam) - 1)
        seen = set(); Rs = []; worstH = 0; worstA = 1e9; worst7 = 0
        rho_ok = True
        for pos in itertools.combinations(range(1, k), m-1):
            w = [0]*k; w[0] = 1
            for i in pos: w[i] = 1
            rots = [tuple(w[i:] + w[:i]) for i in range(k) if w[i] == 1]
            key = min(rots)
            if key in seen: continue
            seen.add(key)
            xs = [Fr(cvalue(list(rw)), D) for rw in rots]
            xmin = min(xs)
            H = len(xs)/sum(1/v for v in xs); A = sum(xs)/len(xs)
            worstH = max(worstH, float(H)/up)
            worstA = min(worstA, float(A) - (up - 1/6))
            rho = Fr(4, 3) - 1/(3*xmin)
            bound = 1 - Fr(m//2, m)*(1 - 1/rho)
            if xmin >= 1:      # Proposition means (iii) assumes x_min >= 1
                if not xmin < float(bound)*up: rho_ok = False
                worst7 = max(worst7, float(xmin)/(float(bound)*up))
            Rs.append(float(xmin)/up)
        Rs = np.array(Rs)
        if m == 11: R11[k] = (Rs, low/up)
        cls_stats.append((m, k, len(seen), Rs.min(), Rs.max(), low/up, worstH, worstA, rho_ok, worst7))
say('  m  k  classes   R=xmin/U range         L/U     max H/U   min A-(U-1/6)  prop(iii) ok (x_min>=1)  max xmin/bound')
for (m, k, n, rmin, rmax, lu, wh, wa, ok, w7) in cls_stats:
    say(f'  {m:2d} {k:2d} {n:7d}  [{rmin:.4f}, {rmax:.4f}]   {lu:.4f}   {wh:.4f}     {wa:.4f}        {ok}     {w7:.4f}')
say('  overall: max H/U = %.4f, min A-(U-1/6) = %.4f, max R = %.4f, prop (iii) holds in all: %s'
    % (max(s[6] for s in cls_stats), min(s[7] for s in cls_stats), max(s[4] for s in cls_stats),
       all(s[8] for s in cls_stats)))
say('  min over (m,k) of (R_min / (L/U)) =', min(s[3]/s[5] for s in cls_stats), ' (1 means the lower end is attained)')
nt = [s for s in cls_stats if s[1] < 2*s[0]]
say('  max R with k < 2m: %.4f' % max(s[4] for s in nt))
top = max(cls_stats, key=lambda s: s[4]); topnt = max(nt, key=lambda s: s[4])
say('  max R attained at (m, k) = (%d, %d); with k < 2m at (m, k) = (%d, %d); 3 log(4/3) = %.4f'
    % (top[0], top[1], topnt[0], topnt[1], 3*math.log(4/3)))

# ---------------- Part D ----------------
say('\nPart D: balanced (Christoffel) words with k = ceil(m theta), relative to U = m/(3 Lambda) = F(m)')

def balanced_cycle(m, k):
    """Least odd, greatest odd and greatest element of the rational cycle of the
    balanced word, as numerators over D = 2^k - 3^m (D is odd, so the parity of
    x = n/D is that of n; T sends n to (3n + D)/2 or n/2)."""
    w = [(i+1)*m//k - i*m//k for i in range(k)]
    D = 2**k - 3**m
    n = cvalue(w); n0 = n; odd = []; top = n
    for t in range(k):
        assert n % 2 == w[t]
        if n % 2:
            odd.append(n); n = (3*n + D)//2
        else:
            n //= 2
        top = max(top, n)
    assert n == n0
    return min(odd), max(odd), top, D

def balanced_ratios(m):
    k = int(mp.ceil(m*TH))
    lo_, hi_, top, D = balanced_cycle(m, k)
    U = m/(3*(k*L2 - m*L3))
    return k, float(mp.mpf(lo_)/D/U), float(mp.mpf(hi_)/D/U), float(mp.mpf(top)/D/U)

say('  the m of the first version of this computation (gcd(m, k) > 1 means a repeated word):')
for m in [5, 10, 20, 50, 100, 200, 400, 800]:
    k, a1, a2, a3 = balanced_ratios(m)
    say(f'   m={m:4d} k={k:4d} gcd={math.gcd(m, k):3d}  xmin/F={a1:.5f}  max odd/F={a2:.5f}  max/F={a3:.5f}')
prim = []
for m in range(100, 2001):
    k = int(mp.ceil(m*TH))
    if math.gcd(m, k) == 1:
        prim.append((m,) + balanced_ratios(m))
for lo_m, hi_m in [(100, 800), (100, 2000), (801, 2000)]:
    sel = [p for p in prim if lo_m <= p[0] <= hi_m]
    say('  primitive words, %d <= m <= %d (%d words): xmin/F in [%.5f, %.5f], max odd/F in [%.5f, %.5f], max/F in [%.5f, %.5f]'
        % (lo_m, hi_m, len(sel), min(p[2] for p in sel), max(p[2] for p in sel), min(p[3] for p in sel),
           max(p[3] for p in sel), min(p[4] for p in sel), max(p[4] for p in sel)))
say('  1/(2 ln 2) =', 1/(2*math.log(2)), ' 1/ln 2 =', 1/math.log(2))

# ---------------- figure ----------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(WIDTH, 2.55),
                               gridspec_kw=dict(width_ratios=[1.25, 1], wspace=0.62))
ks = sorted(R11)
cm = plt.get_cmap('plasma')
norm = BoundaryNorm(np.arange(ks[0]-0.5, ks[-1]+1.5), cm.N)
bins = np.linspace(-2.0, 0.05, 42)
for k in ks:
    Rs, lu = R11[k]
    h, e = np.histogram(np.log10(Rs), bins=bins, density=True)
    ax1.step(e[:-1], h, where='post', color=cm(norm(k)), lw=1.1)
    ax1.plot([math.log10(lu)]*2, [0, 0.35], color=cm(norm(k)), lw=1.0, ls=':')
ax1.axvline(0, color=PAL['ink'], lw=1.0, ls='--')
ax1.set_ylim(0, ax1.get_ylim()[1]*1.12)
ax1.text(-0.04, ax1.get_ylim()[1]*0.95, 'upper end', ha='right', va='top', fontsize=6.5, color=PAL['ink'])
ax1.text(-1.97, 0.45, 'lower ends', fontsize=6.5, color=PAL['plum'])
ax1.set_xlabel(r'$\log_{10}(x_{\min}/U)$, $m=11$')
ax1.set_ylabel('density')
sm = plt.cm.ScalarMappable(norm=norm, cmap=cm); sm.set_array([])
cb = fig.colorbar(sm, ax=ax1, ticks=ks, pad=0.02, aspect=25); cb.set_label(r'$k$')
cb.ax.tick_params(labelsize=6.5)
pts = [(s[0], s[11], s[1] - math.ceil(s[0]*float(TH))) for s in summary
       if not math.isnan(s[11]) and s[3] >= 50]
sc2 = ax2.scatter([p[0] for p in pts], [p[1] for p in pts], c=[p[2] for p in pts],
                  cmap=CMAP_SEQ, s=13, edgecolor=PAL['ink'], linewidth=0.3, zorder=3)
ax2.axhline(1, color=PAL['coral'], lw=1.0, ls='--')
ax2.set_yscale('log'); ax2.set_xlabel(r'number of ones $m$')
ax2.set_ylabel(r'closest approach / uniform model')
cb2 = fig.colorbar(sc2, ax=ax2, pad=0.02, aspect=25); cb2.set_label(r'$k-\lceil m\theta\rceil$')
cb2.ax.tick_params(labelsize=6.5)
fig.savefig(BOOK + '/figures/ch8_words.pdf'); plt.close(fig)
