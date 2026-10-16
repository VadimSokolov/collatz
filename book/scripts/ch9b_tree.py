"""The backward Syracuse tree from 1 below 2^24 for the chapter on positive density from the backward tree (adapted from
N. Polson's tree.py). Every odd n < 2^24 reaches 1, so each is the source of exactly one word valid at the target 1 that
meets 1 only at its end (since S(1) = 1, the words (2, ..., 2, w) are valid too and are not counted): its
letters are the valuations v_2(3y+1) along the Syracuse orbit of n down to 1, its length l(n) is the number of Syracuse
steps, its degree g(n) the sum of the letters, and its weight 3^l(n) / 2^g(n). The orbits are followed forward until they
drop below their start, and the counts of the smaller odd integer are reused. For each window [X, 32X) with X = 2^4, ...,
2^19: the number of odd integers, the number whose orbit stays below 2^24 (the sources that a breadth first search of the
tree below 2^24 finds, as in Polson's table), the total weight beta_X(1) over the window and over the sources the search
finds, and count/(X beta_X). About 2 minutes and 1 GB. Writes data/ch9b_tree.json."""
import os, json, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, 'data'); os.makedirs(DATA, exist_ok=True)
LIMIT = 2**24; H = LIMIT // 2                       # odd n = 2i + 1 with i < H
L = [0] * H; G = [0] * H; low = [True] * H          # n = 1 is the source of the empty word
for i in range(1, H):
    n = 2 * i + 1; y = n; l = g = 0; top = n
    while y >= n:
        y = 3 * y + 1; a = (y & -y).bit_length() - 1; y >>= a; l += 1; g += a
        if y > top: top = y
    j = y >> 1; L[i] = l + L[j]; G[i] = g + G[j]; low[i] = top < LIMIT and low[j]
L = np.array(L, dtype=np.int64); G = np.array(G, dtype=np.int64); low = np.array(low)
wt = np.exp(L * math.log(3) - G * math.log(2))
rows = []
for e in range(4, 20):
    X = 2**e; s = slice(X // 2, 16 * X)                 # the odd n in [X, 32X)
    cnt = 16 * X - X // 2; found = int(low[s].sum()); beta = float(wt[s].sum()); beta_found = float(wt[s][low[s]].sum())
    rows.append(dict(e=e, odd=cnt, found=found, beta=beta, beta_found=beta_found, ratio=cnt / (X * beta),
                     ratio_found=found / (X * beta_found), maxlen=int(L[s].max()), maxwt_times_n=float((wt[s] * (2 * np.arange(X // 2, 16 * X) + 1)).max())))
    print(rows[-1], flush=True)
json.dump(dict(limit=LIMIT, rows=rows, maxlen=int(L.max()), nmaxlen=int(2 * int(L.argmax()) + 1)), open(os.path.join(DATA, 'ch9b_tree.json'), 'w'))
