"""Chapter 5b: growth of the doubling chain construction at large scales.
R_n(X) is the count of ch5b_counts.py (root n, its doubling chain, and recursively the
children (2^(j+1) p - 1)/3 of the chain elements 2^j p with j >= 1 and 2^j p = 2, 8 mod 9),
restricted to integers <= X n.  An element is counted at scale X exactly when every chain
element on its path from the root is <= X n, so one traversal at the largest scale gives
R_n(X) for all X.  Prints R_n(2^i) and least squares slopes of log R against log X over the
points of the last decade up to 2^20 and up to 2^36 (the four points 2^17..2^20 and
2^33..2^36), for the roots n = 4, 5, 7.  Run: python3 -I ch5b_chains.py"""
import numpy as np

IMAX = 36
ROOTS = [4, 5, 7]

def peaks(n, top):
    """path maxima of all elements of the construction from root n with path maximum <= top"""
    out = []; stack = [(n, n)]
    while stack:
        r, pk_r = stack.pop()
        p = r; j = 0
        while p <= top:
            pk = max(pk_r, p)
            out.append(pk)
            if j >= 1 and p % 9 in (2, 8):
                stack.append(((2*p - 1)//3, pk))
            p *= 2; j += 1
    return np.sort(np.array(out, dtype=np.float64))

if __name__ == '__main__':
    Xs = 2.0**np.arange(0, IMAX + 1)
    for n in ROOTS:
        pk = peaks(n, n*2**IMAX)
        R = np.searchsorted(pk, Xs*n, side='right')
        lx, lr = np.log2(Xs), np.log2(R)
        s20 = np.polyfit(lx[17:21], lr[17:21], 1)[0]
        s36 = np.polyfit(lx[33:37], lr[33:37], 1)[0]
        print(f"n={n}: R(2^20)={R[20]} R(2^36)={R[36]}  slope over 2^17..2^20 {s20:.4f}  "
              f"over 2^33..2^36 {s36:.4f}", flush=True)
        print("   R(2^i), i=0..20:", list(R[:21]), flush=True)
