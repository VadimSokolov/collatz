import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
os.makedirs(BOOK + '/scripts/data', exist_ok=True)  # intermediate results that later scripts load
# Trees of slow integers as a function of the box size B and the level m (project P9).
# Tree = coalescence class: integers in {sigma > m} on [2,2^B) whose orbits merge before the
# remaining stopping time falls to r = m - 100.  Key of n = the orbit element with sigma = r.
import numpy as np, math, json
tau = np.load(BOOK + '/scripts/data/ch3_tau_small.npy'); sig = tau.astype(np.int64); M = len(sig)
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
out = []
for B in (18, 20, 22, 24, 26):
    for m in (200, 250, 300, 350, 400):
        S = np.nonzero(sig[2:1 << B] > m)[0] + 2
        if len(S) == 0: out.append((B, m, 0, 0, 0.0, 0)); continue
        keys = np.array([key(int(n), m - 100) for n in S])
        u, c = np.unique(keys, return_counts=True)
        out.append((B, m, int(len(S)), int(len(u)), float(c.mean()), int(c.max())))
        print(f"B={B} m={m}: {len(S):8d} integers, {len(u):6d} trees, mean size {c.mean():8.1f}, largest {c.max():6d}")
json.dump(out, open(BOOK + '/scripts/data/ch10c_trees.json', 'w'))
