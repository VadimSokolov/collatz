"""Exact test of Theorem 20.2 on random rational cycles (fractions, no rounding).
For a word k_0..k_{r-1} with sum s, 2^s > 3^r, the rational cycle has
x_0 = c/(2^s - 3^r). We rotate to the minimum and test
3 M Lambda < B_r(M) = sum_t min(1, 2^{floor(t theta + e_t) - t theta})."""
import random, math, mpmath as mp
from fractions import Fraction
mp.mp.dps = 50
TH = mp.log(3, 2)
def cycle(word):
    r, s = len(word), sum(word)
    K = [0]
    for k in word[:-1]: K.append(K[-1]+k)
    c = sum(3**(r-1-j) * 2**K[j] for j in range(r))
    x = Fraction(c, 2**s - 3**r); xs = []
    for k in word:
        xs.append(x); x = (3*x + 1) / 2**k
    assert x == xs[0]
    return xs
def B(r, M):
    u = mp.log(1 + 1/(3*mp.mpf(M)), 2)
    tot = mp.mpf(0)
    for t in range(r):
        v = t*TH; e = t*u
        tot += min(1, mp.mpf(2)**(mp.floor(v + e + mp.mpf(10)**-30) - v))
    return tot
random.seed(1)
worst = 0; n = 0
for trial in range(3000):
    r = random.randint(1, 25)
    s = math.floor(r*float(TH)) + random.randint(1, 3)
    if 2**s <= 3**r: continue
    # random composition of s into r positive parts
    cuts = sorted(random.sample(range(1, s), r-1)) if r > 1 else []
    w = [b-a for a, b in zip([0]+cuts, cuts+[s])]
    xs = cycle(w); i = min(range(r), key=lambda j: xs[j])
    w = w[i:]+w[:i]; M = xs[i]
    Lam = s*mp.log(2) - r*mp.log(3)
    lhs = 3*mp.mpf(M.numerator)/M.denominator*Lam
    rhs = B(r, float(M))
    worst = max(worst, float(lhs/rhs)); n += 1
print("cycles tested:", n, " max of 3 M Lambda / B_r(M):", worst)

random.seed(1)
for trial in range(3000):
    r = random.randint(1, 25)
    s = math.floor(r*float(TH)) + random.randint(1, 3)
    if 2**s <= 3**r: continue
    cuts = sorted(random.sample(range(1, s), r-1)) if r > 1 else []
    w = [b-a for a, b in zip([0]+cuts, cuts+[s])]
    xs = cycle(w); i = min(range(r), key=lambda j: xs[j])
    M = xs[i]; Lam = s*mp.log(2) - r*mp.log(3)
    lhs = 3*mp.mpf(M.numerator)/M.denominator*Lam; rhs = B(r, float(M))
    if lhs > rhs:
        print(r, s, float(M), float(Lam), float(lhs), float(rhs), w[i:]+w[:i], [float(v) for v in xs[i:]+xs[:i]][:6]); break
