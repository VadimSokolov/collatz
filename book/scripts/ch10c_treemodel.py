import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
os.makedirs(BOOK + '/scripts/data', exist_ok=True)  # intermediate results that later scripts load
# Trees at level m on [2,2^B): roots, the model prediction for their sizes, and the delay records.
# Model: a backward path of length d with j odd steps from the root x exists with probability 3^-j
# (integrality), and the predecessor it reaches has size about x 2^d / 3^j; it belongs to the tree
# if d > 100 (so that sigma > m) and the size is below 2^B.
import numpy as np, math, json
from math import comb, log2
B = 26; M = 1 << B
tau = np.load(BOOK + '/scripts/data/ch3_tau_small.npy'); sig = tau.astype(np.int64)
def T(x): return (3 * x + 1) // 2 if x & 1 else x // 2
def sigma(x):
    k = 0
    while True:
        if x < M: return sig[x] + k
        x = T(x); k += 1
def key(n, r):
    # the orbit point with remaining stopping time r, followed above the box when necessary
    x = n
    while True:
        s = sigma(x)
        if s == r: return x
        if s < r: return -1
        x = T(x)
def predicted(x, dmin=100, dmax=700):
    tot = 0.0; lx = log2(x)
    for d in range(dmin + 1, dmax):
        for j in range(d + 1):
            if lx + d - j * log2(3) < B: tot += comb(d, j) * 3.0 ** (-j)
    return tot
# delay records (under T) below 2^26
rec = []; best = -1
for n in range(2, M):
    if sig[n] > best: best = sig[n]; rec.append((n, int(best)))
recset = {n for n, _ in rec}
out = {}
for m in (350, 400, 440):
    r = m - 100
    S = np.nonzero(sig[2:] > m)[0] + 2
    keys = np.array([key(int(n), r) for n in S])
    rows = []
    for x in np.unique(keys):
        members = S[keys == x]
        pred = predicted(int(x)) if x > 0 else float('nan')
        top = members[np.argmax(sig[members])]
        rows.append((int(x), int(len(members)), pred, int(top), int(sig[top]), bool(len(set(members.tolist()) & recset) > 0)))
    rows.sort(key=lambda t: -t[1])
    out[m] = rows
    preds = np.array([t[2] for t in rows]); sizes = np.array([t[1] for t in rows])
    print(f"m={m}: {len(rows)} trees; model prediction max {preds.max():.3f}, median {np.median(preds):.3g}; "
          f"observed size / prediction: median {np.median(sizes / preds):.3g}, min {np.min(sizes / preds):.3g}")
    for x, size, pred, top, st, isrec in rows[:10]:
        print(f"   root {x:12d} (log2 {(log2(x) if x > 0 else float('nan')):5.2f})  size {size:5d}  model {pred:9.3g}  deepest member {top:9d} sigma {st}  contains delay record: {isrec}")
json.dump({str(k): v for k, v in out.items()}, open(BOOK + '/scripts/data/ch10c_treemodel.json', 'w'))
print("delay records below 2^26 with sigma>350:", [(n, s) for n, s in rec if s > 350])
