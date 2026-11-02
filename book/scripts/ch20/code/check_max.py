"""Exact test of Theorem 20.5 (maximum element) on random rational cycles:
3 N Lambda >= (1 + 1/(3M))^{-1} e^{-Lambda} (2 B_r - 1),  B_r = sum_{t<r} 2^{-{t theta}}."""
import random, math, mpmath as mp
from check_random import cycle, TH
mp.mp.dps = 50
random.seed(2); worst = 9; n = 0
for trial in range(3000):
    r = random.randint(1, 25)
    s = math.floor(r*float(TH)) + random.randint(1, 3)
    if 2**s <= 3**r: continue
    cuts = sorted(random.sample(range(1, s), r-1)) if r > 1 else []
    w = [b-a for a, b in zip([0]+cuts, cuts+[s])]
    xs = cycle(w); N = max(xs); M = min(xs)
    Lam = s*mp.log(2) - r*mp.log(3)
    B = sum(mp.mpf(2)**(-mp.frac(t*TH)) for t in range(r))
    lhs = 3*mp.mpf(N.numerator)/N.denominator*Lam
    rhs = 1/(1 + 1/(3*mp.mpf(M.numerator)/M.denominator))*mp.e**(-Lam)*(2*B - 1)
    worst = min(worst, float(lhs/rhs)); n += 1
print("cycles:", n, " min of lhs/rhs:", worst)
