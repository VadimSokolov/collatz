"""Numbers of chapter 20 added by the review of 10 November 2026.
(1) Orbits of B_m(x) = x/3, (x-1)/3, (m x - 2)/3 by x mod 3, from 1 <= n <= 20000: how many enter 0, a nontrivial cycle, or pass 10^60.
(2) The sharp bound B_17/(3(2^delta-1)) at the rung r = 17, s = 27.
(3) The largest minimum of a rational cycle against the bound B'_r(delta_c)/(3(2^delta_c - 1)) for r = 3, 6, 9, 12 with slope 5/3, by enumeration.
(4) The share of starts landing in cycles with 85 odd steps, to two decimals, from basin_17_27.pkl (written by basin.py 17 27 160000 80), if present.
Usage: python review_checks.py"""
import math, itertools, os, pickle, collections
th = math.log2(3)

def orbit_class(m, n, cap=10**60):
    seen = set(); x = n
    while x not in seen:
        if x == 0: return 'zero'
        if x > cap: return 'escaped'
        seen.add(x); r = x % 3
        x = x // 3 if r == 0 else (x - 1) // 3 if r == 1 else (m * x - 2) // 3
    return 'zero' if x == 0 else 'cycle'

for m in (7, 13, 19):
    c = collections.Counter(orbit_class(m, n) for n in range(1, 20001))
    print('B_%d:' % m, dict(c))

r, s = 17, 27
d = s - r * th
Br = sum(2 ** (-((t * th) % 1)) for t in range(r))
print('r=17: B_r/(3(2^delta-1)) = %.4f' % (Br / (3 * (2 ** d - 1))))

def phi(w): return sum(w[j] * 3 ** sum(w[j + 1:]) * 2 ** j for j in range(len(w)))
def Bp(r, dc): return sum(2 ** (math.floor(t * th + dc) - t * th) for t in range(r))
for r, s in ((3, 5), (6, 10), (9, 15), (12, 20)):
    best = 0
    for pos in itertools.combinations(range(s), r):
        w = [0] * s
        for p in pos: w[p] = 1
        best = max(best, min(phi(w[k:] + w[:k]) for k in range(s)))
    dc = s - r * th
    print('(r,s)=(%d,%d): largest minimum %.4f, bound %.4f' % (r, s, best / (2 ** s - 3 ** r), Bp(r, dc) / (3 * (2 ** dc - 1))))

if os.path.exists('basin_17_27.pkl'):
    dd, rows = pickle.load(open('basin_17_27.pkl', 'rb'))
    tot = sum(x[1] for x in rows)
    by = collections.Counter()
    for x in rows: by[x[2]] += x[1]
    for k in sorted(by)[:6]: print('odd steps %d: share %.2f%%' % (k, 100 * by[k] / tot))
