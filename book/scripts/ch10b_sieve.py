import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
os.makedirs(BOOK + '/scripts/data', exist_ok=True)  # intermediate results that later scripts load
# The residue sieve: how many classes mod 2^k are not certified to descend within k steps of T
# (coefficient stopping time > k), computed exactly by the dynamic programme of Chapter 2:
# N_k(m) = (N_{k-1}(m) + N_{k-1}(m-1)) [3^m > 2^k], alive(k) = sum_m N_k(m).
import math, json
K = 400
N = {0: 1}; rows = []
for k in range(1, K + 1):
    new = {}
    for m, c in N.items():
        for mm in (m, m + 1):
            if 3 ** mm > 2 ** k:
                new[mm] = new.get(mm, 0) + c
    N = new
    alive = sum(N.values())
    rows.append((k, alive))
out = []
for k, a in rows:
    frac = a / 2 ** k
    out.append((k, a, frac, math.log2(a) / k if a > 0 else 0.0))
for k in (10, 20, 40, 60, 80, 100, 150, 200, 300, 400):
    k_, a, f, e = out[k - 1]
    print(f"k={k:4d} uncovered classes = {a:.4e}  proportion = {f:.4e}  log2(count)/k = {e:.4f}")
# Chernoff rate D: proportion ~ exp(-D k) up to polynomial factors
k1, k2 = 200, 400
D = -(math.log(out[k2 - 1][2]) - math.log(out[k1 - 1][2])) / (k2 - k1)
print("empirical decay rate of proportion between k=200 and 400:", D, " growth exponent of count 1 - D/ln2 =", 1 - D / math.log(2))
json.dump(out, open(BOOK + '/scripts/data/ch10b_sieve.json', 'w'))
