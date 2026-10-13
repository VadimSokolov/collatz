import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
os.makedirs(BOOK + '/scripts/data', exist_ok=True)  # intermediate results that later scripts load
# MM approach pressed harder: (1) residue graph cycle means and positive cycles, (2) block MM over runs, (3) run length potential.
import numpy as np, math, json
from fractions import Fraction
L32, L2, L3 = math.log(1.5), math.log(2), math.log(3)
def T(n): return n // 2 if n % 2 == 0 else (3 * n + 1) // 2
# ---------- (1) residue graph G_d mod 2^d: nodes a, edges a -> lifts of T(a) mod 2^(d-1) to mod 2^d, weight log(3/2) if a odd else -log 2.
# Maximum cycle mean by Karp's algorithm.
def karp(d):
    M = 2 ** d; W = {}
    for r in range(M):
        t = T(r) % (M // 2)
        for lift in (t, t + M // 2):
            W.setdefault(r, []).append((lift, L32 if r % 2 else -L2))
    # Karp: D[m][v] = max weight of walk of length m ending at v from any start (use a super source)
    INF = -1e18
    D = np.full((M + 1, M), INF); D[0, :] = 0.0
    for m in range(1, M + 1):
        for u in range(M):
            if D[m - 1, u] > INF / 2:
                for v, w in W[u]:
                    if D[m - 1, u] + w > D[m, v]: D[m, v] = D[m - 1, u] + w
    best = -1e18
    for v in range(M):
        if D[M, v] > INF / 2:
            val = min((D[M, v] - D[m, v]) / (M - m) for m in range(M) if D[m, v] > INF / 2)
            best = max(best, val)
    return best
for d in range(2, 11):
    print(f"d={d}: max cycle mean = {karp(d):.6f}   (log(3/2) = {L32:.6f})")
# all simple cycles of G_d, each found once from its least residue by depth first search
def wt(a): return L32 if a % 2 else -L2
def simple_cycles(d):
    M = 2 ** d
    succ = {a: sorted({T(a) % (M // 2), T(a) % (M // 2) + M // 2}) for a in range(M)}
    out = []
    def dfs(s, v, path, onpath):
        for u in succ[v]:
            if u == s: out.append(list(path))
            elif u > s and u not in onpath:
                onpath.add(u); path.append(u); dfs(s, u, path, onpath); path.pop(); onpath.discard(u)
    for s in range(M): dfs(s, s, [s], {s})
    return out
for d in (3, 4, 5):
    cyc = simple_cycles(d)
    pos = [c for c in cyc if sum(wt(a) for a in c) > 0]
    print(f"G_{d}: {len(cyc)} simple cycles, {len(pos)} of positive weight")
cyc = simple_cycles(5)
pos = sorted((c for c in cyc if sum(wt(a) for a in c) > 0), key=lambda c: -sum(wt(a) for a in c) / len(c))
for c in pos[:4]: print("  G_5 positive cycle, mean %.4f, residues %s" % (sum(wt(a) for a in c) / len(c), sorted(c)))
print("largest cycle mean of G_5 without the loop at -1: %.6f; log(243/64)/6 = %.6f"
      % (max(sum(wt(a) for a in c) / len(c) for c in cyc if c != [31]), math.log(243 / 64) / 6))
# the rational cycle of T through -19/11 (parity word 1110) traces the positive 4 cycle {7, 9, 14, 27} of G_5
x = Fraction(-19, 11); orb = []
for _ in range(4):
    orb.append(x); x = x / 2 if x.numerator % 2 == 0 else (3 * x + 1) / 2
print("cycle of -19/11:", [str(z) for z in orb], "back to start:", x == orb[0], "residues mod 32:", sorted(z.numerator * pow(z.denominator, -1, 32) % 32 for z in orb))
# the negative cycles as residue cycles: weights
for cyc in ([-1], [-5, -7, -10], [-17, -25, -37, -55, -82, -41, -61, -91, -136, -68, -34]):
    r = sum(1 for x in cyc if x % 2); s = len(cyc)
    print("negative cycle", cyc[0], "odd", r, "steps", s, "total weight", r * L3 - s * L2, "mean", (r * L3 - s * L2) / s,
          "residues mod 32", [x % 32 for x in cyc])
# ---------- (2) block MM: x = 2^j q - 1 (q odd), run of j-1 ascents then one descent; block image (3^j q - 1)/2^{v}, v = v2(3^j q - 1)
def v2(n):
    c = 0
    while n % 2 == 0: n //= 2; c += 1
    return c
def S(x):
    y = 3 * x + 1
    while y % 2 == 0: y //= 2
    return y
# exact check of the block formula, of the identity q(2^(j+v) - 3^j) = 2^v (2^j q - u) - 1 with u = S^(j)(x), and of the
# descent criterion v > j log2(3/2), over odd q < 2^16 for j = 1..12; the only block fixed point is x = 1 (Steiner)
res = {}
for j in range(1, 13):
    cnt = 0; desc = 0; tot = 0; cond = {}
    for q in range(1, 2 ** 16, 2):
        x = 2 ** j * q - 1
        y = x
        for _ in range(j): y = S(y)
        v = v2(3 ** j * q - 1)
        assert y == (3 ** j * q - 1) // 2 ** v, (j, q)
        assert q * (2 ** (j + v) - 3 ** j) == 2 ** v * (2 ** j * q - y) - 1, (j, q)
        cnt += 1; desc += (y < x)
        need = v > j * math.log2(1.5)
        assert (y < x) == need or x == 1, (j, q, y, x, v)
        if y == x: print("block fixed point: x =", x, "j =", j, "v =", v)
    thr = math.floor(j * math.log2(1.5)) + 1          # smallest v that descends
    res[j] = (desc / cnt, 2.0 ** (-(thr - 1)))
    print(f"j={j:2d}: block descends with frequency {desc/cnt:.5f}; model 2^-(ceil(j log2 1.5)-1) = {2.0**(-(thr-1)):.5f}; threshold v >= {thr}")
json.dump(res, open(BOOK + '/scripts/data/ch5c_blocks.json', 'w'))
# ---------- (3) run length potential V_c(x) = log x + c v2(x+1): fraction of odd 3 <= x < N with V_c(S(x)) >= V_c(x)
# An ascending stride from x fails exactly when log(3/2 + 1/(2x)) >= c, so for c > log(3/2) only finitely many fail,
# and none when c > log(5/3) (the stride 3 -> 5); the failures are listed explicitly.
N = 2 ** 20
data = []
for x in range(3, N, 2):
    y = S(x); data.append((x, y, math.log(y / x), v2(x + 1), v2(y + 1)))
for c in (0.45, 0.5, math.log(5 / 3), 0.52, 1.0, 2.0):
    fail = 0; asc = []
    for x, y, lr, a, b in data:
        if lr + c * (b - a) >= 0:
            fail += 1
            if y > x: asc.append(x)
    print(f"c={c:.4f}: potential fails at fraction {fail/len(data):.4f} of odd x; ascending failures at x = {asc}")
# model: a descending stride with k halvings fails iff it lands on y with v2(y+1) - 1 >= (k log 2 - log 3 - eps(x))/c,
# v2(y+1) ~ Geometric(1/2) on {1, 2, ...} independent of k: P(v2(y+1) - 1 >= t) = 2^-t
def pmodel(c): return sum(2.0 ** (-k) * 2.0 ** (-math.ceil((k * L2 - L3) / c)) for k in range(2, 80))
for c in (1.0, 2.0, L32 + 1e-9):
    print(f"model failure fraction for c={c:.6f}: {pmodel(c):.6f}")
