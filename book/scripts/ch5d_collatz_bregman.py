# Bregman (Itakura and Saito) decomposition of the Syracuse step: for M = S(x)/x,
#   log M = (M - 1) - D(M),  D(t) = t - 1 - log t >= 0.
# (1) uniform dissipation constant c* = inf over odd x >= 3 of D(S(x)/x), and its large x limit;
# (2) the model: E[M] = 1 and E[D] = log(4/3) for M = 3/2^k, k ~ Geometric(1/2);
# (3) orbit statistics from random large odd starting values: mean D per step, mean (M-1), and
#     the variance of M - 1 against the model value 2/7.
import math, random, json, os
D = lambda t: t - 1 - math.log(t)
def S(x):
    y = 3 * x + 1
    while y % 2 == 0: y //= 2
    return y
best = (9e9, None)
for x in range(3, 2_000_001, 2):
    d = D(S(x) / x)
    if d < best[0]: best = (d, x)
print("c* =", best[0], "attained at x =", best[1], " exact D(7/9) =", D(7 / 9))
print("large x limit D(3/4) =", D(0.75), "; D(3/2) =", D(1.5), "; D(3/8) =", D(3 / 8))
# minimum of D over x >= X for a few X
for X in (10, 100, 10 ** 4, 10 ** 6):
    m = min(D(S(x) / x) for x in range(X | 1, X + 400001, 2)); print(f"  inf over x in [{X}, {X}+4e5]: {m:.6f}")
EM = sum(3 / 2 ** k * 2 ** -k for k in range(1, 80)); ED = sum(D(3 / 2 ** k) * 2 ** -k for k in range(1, 80))
VM = sum((3 / 2 ** k - 1) ** 2 * 2 ** -k for k in range(1, 80))
print("model: E M =", EM, " E D =", ED, " log(4/3) =", math.log(4 / 3), " Var M =", VM, " 2/7 =", 2 / 7)
random.seed(7); stats = []
for trial in range(400):
    x = random.getrandbits(3000) | 1; x0 = x; n = 0; sD = sL = sL2 = 0.0; path = []
    while x > 1 and n < 40000:
        y = S(x); M = y / x if x < 10 ** 300 else math.exp(math.log(y) - math.log(x)) if False else None
        # exact ratio via logs of big integers
        lm = (y.bit_length() - x.bit_length()) * math.log(2) + math.log((y >> max(0, y.bit_length() - 60)) / (1 << min(60, y.bit_length()))) - math.log((x >> max(0, x.bit_length() - 60)) / (1 << min(60, x.bit_length())))
        M = math.exp(lm); d = M - 1 - lm
        sD += d; sL += M - 1; sL2 += (M - 1) ** 2; n += 1; x = y
    stats.append((n, sD / n, sL / n, sL2 / n - (sL / n) ** 2))
n_tot = sum(s[0] for s in stats)
mD = sum(s[0] * s[1] for s in stats) / n_tot; mL = sum(s[0] * s[2] for s in stats) / n_tot; vL = sum(s[0] * s[3] for s in stats) / n_tot
print(f"orbits from 3000 bit starts: {len(stats)} orbits, {n_tot} Syracuse steps; mean D = {mD:.5f} (model {math.log(4/3):.5f}); mean M-1 = {mL:.5f} (model 0); var M = {vL:.5f} (model {2/7:.5f})")
json.dump(dict(cstar=best[0], cstar_x=best[1], D34=D(0.75), ED=ED, EM=EM, VM=VM, orbit_meanD=mD, orbit_meanL=mL, orbit_varL=vL, steps=n_tot),
          open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'ch5d_bregman.json'), 'w'))
