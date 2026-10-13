"""Chapter 5b: the certificates of ch5b_tree.py checked in exact integer arithmetic.
For each exponent g of the table of thm:ch5b:exponent, and for g = 0.47705 at k = 7 (cor:ch5b:count), the vector w
found by power iteration in ch5b_tree.py is rounded to positive integers, and a lower bound for
min_r (Phi_g w)_r / w_r is computed with integers only: 3^g and 2^(-g) are rounded down at 2^-S, every product is
rounded down, the sum over j is taken over the first period j = 1..P, and the factor 1/(1 - 2^(-Pg)) that the later periods
contribute is replaced by 1/(1 - q_P) with q_P <= 2^(-Pg).  Each of these can only lower the bound, so a printed bound
above 1 certifies the exponent.
Run: python3 -I ch5b_tree_exact.py"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction
from mpmath import mp, mpf, floor as mfloor
from ch5b_tree import build, certify

S = 256                                    # fixed point: a real x >= 0 is held as an integer X <= x 2^S
mp.prec = 600

def lower(x):
    return int(mfloor(x * mpf(2)**S)) - 1  # the -1 absorbs the error of the 600 bit evaluation

def margin(k, g):
    classes, lifts, P = build(k); n = len(classes)
    _, c = certify(lifts, n, P, float(g), 2000)
    w = [int(round(v * 2**60)) for v in c]
    assert min(w) > 0
    gm = mpf(g.numerator) / g.denominator
    q = lower(mpf(2)**(-gm)); t3 = lower(mpf(3)**gm)
    qp = [0, q]                            # qp[i] <= 2^(-g i) 2^S
    for i in range(2, P + 2):
        qp.append(qp[-1] * q >> S)
    best = None
    for r in range(n):
        lo = None
        for rows, cols, jj in lifts:       # one lift of r: the sum over its fertile j
            sel = rows == r
            s = sum(qp[int(j) + 1] * w[int(col)] for col, j in zip(cols[sel], jj[sel]))
            lo = s if lo is None else min(lo, s)
        ratio = Fraction(t3 * lo, w[r] << (2 * S))
        best = ratio if best is None else min(best, ratio)
    return n, best * Fraction(1 << S, (1 << S) - qp[P])   # the later periods: 1/(1 - q_P) <= 1/(1 - 2^(-Pg))

if __name__ == '__main__':
    cases = [(2, '0.416'), (3, '0.435'), (4, '0.449'), (5, '0.460'), (6, '0.470'), (7, '0.477'), (7, '0.47705')]
    for k, g in cases:
        n, m = margin(k, Fraction(g))
        lo9 = Fraction(int(m * 10**9), 10**9)     # rounded down
        print(f"k={k} classes={n} g={g}: min (Phi_g w)_r / w_r >= {float(lo9):.9f}  certified: {m > 1}", flush=True)
