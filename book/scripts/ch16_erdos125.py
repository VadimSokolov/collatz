"""P7: Erdos 125 beyond 4^13.  A = integers with base 3 digits in {0,1}, B = integers with base 4 digits in
{0,1}.  Exact counting function of A+B on [0, N] for N = 4^E using a bit array (N/8 bytes) and compiled
loops.  Reports the density at powers of 4 and 3, the density just below each power of 3 (where the
minima sit), the quality of the nearest approximation 3^m ~ 4^k, and the minimum over [10^4, N].
Usage: python3 ch16_erdos125.py E"""
import numpy as np, numba as nb, sys, json, math, os, time
E = int(sys.argv[1]) if len(sys.argv) > 1 else 15
N = 4 ** E
def digitset(base, limit):
    s = np.array([0], dtype=np.int64); p = 1
    while p <= limit:
        s = np.concatenate([s, s + p]); p *= base
    return np.sort(s[s <= limit])
@nb.njit(cache=True)
def mark(A, B, N, bits):
    for bi in range(B.shape[0]):
        b = B[bi]
        for ai in range(A.shape[0]):
            v = A[ai] + b
            if v > N: break
            bits[v >> 3] |= np.uint8(1 << (v & 7))
@nb.njit(cache=True)
def scan(bits, N, queries, lo):
    # counts at sorted query points, and the minimum density over [lo, N]
    out = np.zeros(queries.shape[0], np.int64); qi = 0; c = 0
    best = 2.0; arg = -1
    for v in range(N + 1):
        if (bits[v >> 3] >> (v & 7)) & 1: c += 1
        while qi < queries.shape[0] and queries[qi] == v:
            out[qi] = c; qi += 1
        if v >= lo:
            d = c / (v + 1.0)
            if d < best: best = d; arg = v
    return out, best, arg
t0 = time.time()
A = digitset(3, N); B = digitset(4, N)
bits = np.zeros(N // 8 + 1, dtype=np.uint8)
mark(A, B, N, bits)
print(f"N = 4^{E} = {N}: |A| = {len(A)}, |B| = {len(B)}, marked in {time.time()-t0:.0f}s", flush=True)
qs = sorted(set([4 ** k for k in range(7, E + 1)] + [3 ** m for m in range(10, 40) if 3 ** m <= N] +
                [3 ** m - 1 for m in range(10, 40) if 3 ** m - 1 <= N]))
cnt, best, arg = scan(bits, N, np.array(qs, dtype=np.int64), 10 ** 4)
dens = {q: c / (q + 1) for q, c in zip(qs, cnt)}
print(f"scanned in {time.time()-t0:.0f}s; minimum density over [1e4, N]: {best:.6f} at N = {arg}", flush=True)
rows = []
for k in range(7, E + 1): print(f"  4^{k}: {dens[4**k]:.4f}")
L43 = math.log(3) / math.log(4)
for m in range(10, 40):
    if 3 ** m - 1 > N: break
    k = round(m * L43); ratio = 3 ** m / 4 ** k
    rows.append(dict(m=m, k=k, ratio=ratio, below=dens[3 ** m - 1], at=dens[3 ** m]))
    print(f"  3^{m}: density just below {dens[3**m-1]:.6f}; nearest 4^{k}, ratio 3^m/4^k = {ratio:.4f}")
json.dump(dict(E=E, N=N, min=best, argmin=int(arg), pow4={str(k): dens[4 ** k] for k in range(7, E + 1)}, pow3=rows),
          open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', f'ch16_erdos125_{E}.json'), 'w'))
