import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
os.makedirs(BOOK + '/scripts/data', exist_ok=True)  # intermediate results that later scripts load
# Trees versus delay records: for each tree at level m, how many members lie on the orbit of its
# deepest member (the trunk) and how many are tributaries; which trees contain a delay record.
import numpy as np, json
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
rec = []; best = -1
for n in range(2, M):
    if sig[n] > best: best = sig[n]; rec.append(n)
recset = set(rec)
res = {}
for m in (300, 350, 400, 440):
    r = m - 100; S = np.nonzero(sig[2:] > m)[0] + 2
    keys = np.array([key(int(n), r) for n in S]); rows = []
    for x in np.unique(keys):
        mem = S[keys == x]; top = int(mem[np.argmax(sig[mem])])
        orbit = set(); y = top
        while sigma(y) > m: orbit.add(y); y = T(y)
        trunk = sum(int(v) in orbit for v in mem)
        rows.append(dict(root=int(x), rootbits=(float(np.log2(x)) if x > 0 else None), size=int(len(mem)), trunk=int(trunk),
                         top=top, topsigma=int(sig[top]), record=bool(set(mem.tolist()) & recset)))
    rows.sort(key=lambda d: -d['size']); res[m] = rows
    nrec = sum(d['record'] for d in rows); big = rows[0]; ntop = sum(d['top'] in recset for d in rows)
    print(f"m={m}: {len(rows)} trees, {nrec} contain a delay record (of sigma under T), in {ntop} the deepest member is one; largest: size {big['size']}, trunk {big['trunk']}, deepest {big['top']} sigma {big['topsigma']}, root {big['root']}")
    for d in rows[:6]: print("   ", d)
json.dump({str(k): v for k, v in res.items()}, open(BOOK + '/scripts/data/ch10c_treerecords.json', 'w'))
