"""Chapter 5b: exact counts pi_n(X) = #{m <= X n : orbit of m under T passes through n}
for several roots n, and counts R_n(X) of the restricted doubling chain construction
(root, doubling chain, children (2^(j+1) n - 1)/3 for j >= 1 when n 2^j = 2,8 mod 9).
Accelerated map T(m) = m/2 (even), (3m+1)/2 (odd)."""
import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
import numpy as np, json

B = 2**25
ROOTS = [4, 5, 7, 8, 11, 13, 16, 17, 19, 20, 22, 23, 25, 26, 28, 29, 31, 32, 34, 35]

def first_below(lo, hi):
    """for m in [lo,hi) return the first iterate of T that is < lo (orbit of m)"""
    m = np.arange(lo, hi, dtype=np.int64); v = m.copy()
    act = np.ones(len(m), bool)
    while act.any():
        w = v[act]
        w = np.where(w % 2 == 0, w//2, (3*w+1)//2)
        v[act] = w
        act[act] = w >= lo
    return v

def labels():
    lab = {n: np.zeros(B+1, bool) for n in ROOTS}
    # small values by direct iteration (all orbits reach 1; stop at cycle {1,2})
    small = 64
    for m in range(1, small):
        seen = []; v = m
        while v not in (1, 2) or len(seen) == 0:
            seen.append(v)
            if v in (1, 2) and len(seen) > 1: break
            v = v//2 if v % 2 == 0 else (3*v+1)//2
            if v in (1, 2): seen.append(v); break
        for n in ROOTS:
            lab[n][m] = n in seen
    lo = small; blk = 2**20
    while lo <= B:
        hi = min(B+1, max(2*lo, lo+1)) if lo < blk else min(B+1, lo+blk)
        F = first_below(lo, hi)
        for n in ROOTS:
            lab[n][lo:hi] = lab[n][F] | (np.arange(lo, hi) == n)
        lo = hi
    return lab

def restricted(n, bound):
    cnt = 0; stack = [n]
    while stack:
        r = stack.pop()
        p = r; j = 0
        while p <= bound:
            cnt += 1
            if j >= 1 and p % 9 in (2, 8):
                c = (2*p - 1)//3
                if c <= bound: stack.append(c)
            p *= 2; j += 1
    return cnt

if __name__ == '__main__':
    lab = labels()
    Xs = [2**i for i in range(0, 21)]
    out = {}
    for n in ROOTS:
        cs = np.cumsum(lab[n])
        N = [int(cs[min(B, X*n)]) for X in Xs if X*n <= B]
        R = [restricted(n, X*n) for X in Xs if X*n <= min(B, 2**31)]
        out[n] = dict(N=N, R=R)
        print(n, n % 9, 'N:', N[-4:], 'R:', R[-4:], flush=True)
    json.dump(dict(Xs=Xs, out=out), open(BOOK + '/scripts/ch5b_counts.json', 'w'))
