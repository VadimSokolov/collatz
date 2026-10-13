import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
os.makedirs(BOOK + '/scripts/data', exist_ok=True)  # intermediate results that later scripts load
# Cluster structure of deep level sets: integers with sigma > m, keyed by the orbit element at which
# the remaining stopping time equals r = m - 100 (orbits that meet at a point of remaining stopping time
# at least r share the key).  Bit flips are counted as ordered pairs (proposals), so each unordered pair twice.
import numpy as np, math, json
B = 26; M = 1 << B
tau = np.load(BOOK + '/scripts/data/ch3_tau_small.npy'); sig = tau.astype(np.int64)
def T(x): return (3 * x + 1) // 2 if x & 1 else x // 2
def sigma(x):
    k = 0
    while True:
        if x < M: return sig[x] + k
        x = T(x); k += 1
def key(n, r):
    x = n
    while True:
        s = sigma(x)
        if s == r: return x
        if s < r: return -1
        x = T(x)
out = {}
for m, r in ((300, 200), (350, 250), (400, 300), (440, 340)):
    S = np.nonzero(sig[2:] > m)[0] + 2
    keys = [key(int(n), r) for n in S]
    uniq, counts = np.unique(keys, return_counts=True)
    # also: fraction of single bit flip neighbours inside the set that share the key
    same = tot = 0
    for n in S:
        for i in range(B):
            n2 = int(n) ^ (1 << i)
            if n2 >= 2 and sig[n2] > m:
                tot += 1; same += (key(n2, r) == key(int(n), r))
    print(f"m={m}: {len(S)} integers, {len(uniq)} clusters (key at remaining time {r}); largest {counts.max()}, "
          f"sizes {sorted(counts.tolist(), reverse=True)[:8]}; in set bit flips (ordered pairs) {tot}, same cluster {same}; "
          f"unordered pairs {tot // 2}, same cluster {same // 2}; crossing flips {tot - same}, per proposal {(tot - same) / (B * len(S)):.2e}")
    out[m] = {'n': int(len(S)), 'clusters': int(len(uniq)), 'sizes': sorted(counts.tolist(), reverse=True), 'nbr': tot, 'same': same}
json.dump(out, open(BOOK + '/scripts/data/ch10c_clusters.json', 'w'))
