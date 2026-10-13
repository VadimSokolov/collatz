"""Chapter 8, Project P1: an exact integer certificate for the sharp cycle bound
r >= 72 057 431 991 at the verified range X = 2^68.  No floating point is used
anywhere in the certificate; every step is an inequality between integers.

Step 1. Rational enclosure lo < theta = log_2 3 < hi, with lo = a/b and
        hi = a'/b', certified by 2^a < 3^b and 3^b' < 2^a'.  The exponents are
        of order 10^11 to 10^12, so the powers are not formed; instead 3^b is
        enclosed between L*2^e and U*2^e by square and multiply with directed
        (floor, ceiling) truncation of the integer mantissa to NBITS bits.
        Each truncation is an integer inequality, so the enclosure is exact.
Step 2. A rational lower bound l2 < ln 2 from the series sum 1/(n 2^n)
        (partial sums are lower bounds), and an upper bound u2 > ln 2 from the
        partial sum plus the tail bound 1/((N+1) 2^N).
Step 3. The partial quotients a_0..a_K of theta: the common prefix of the
        continued fractions of lo and hi (exact Fractions), less the last term.
Step 4. The upper ladder rungs u = q_{k-1} + j q_k (k even), their numerators
        P = p_{k-1} + j p_k, and for every rung u < r_* the exact check
        (P - u*hi) * 3 * X * l2 >= u, i.e. delta(u) >= u/(3 X ln 2).
        For r_* the opposite check (P - r_* lo) * 3 * X * u2 < r_* shows the
        bound cannot be improved by this method; the margin printed for r_* is
        therefore an upper bound, from the other end of the enclosure.
Step 5. Consistency: the delta enclosures of consecutive rungs are strictly
        decreasing, and an exact scan of all r up to SCAN agrees with the ladder.
"""
from fractions import Fraction as Fr
import time, sys

NBITS = 320
X = 2**68

def trunc(m, e, up):
    """Truncate integer mantissa m (>0) to NBITS bits, rounding up or down."""
    sh = m.bit_length() - NBITS
    if sh <= 0:
        return m, e
    if up:
        return -((-m) >> sh), e + sh
    return m >> sh, e + sh

def pow3_enclosure(b):
    """Return (L, eL, U, eU) with L*2^eL <= 3^b <= U*2^eU, integer arithmetic only.
    Square and multiply; after every product the mantissa is truncated to NBITS
    bits, downwards for L and upwards for U, so the inequalities persist."""
    L, eL, U, eU = 1, 0, 1, 0
    for bit in bin(b)[2:]:
        L, eL = trunc(L*L, 2*eL, False); U, eU = trunc(U*U, 2*eU, True)
        if bit == '1':
            L, eL = trunc(3*L, eL, False); U, eU = trunc(3*U, eU, True)
    return L, eL, U, eU

def less_pow2_pow3(a, b):
    """Certify 2^a < 3^b via 2^a < L 2^eL."""
    L, eL, U, eU = pow3_enclosure(b)
    return (a - eL >= 0 and (1 << (a - eL)) < L), (L, eL, U, eU)

def less_pow3_pow2(b, a):
    """Certify 3^b < 2^a via U 2^eU < 2^a."""
    L, eL, U, eU = pow3_enclosure(b)
    return (a - eU >= 0 and U < (1 << (a - eU))), (L, eL, U, eU)

def cf_rational(x, n=200):
    out = []
    for _ in range(n):
        a = x.numerator // x.denominator
        out.append(a)
        f = x - a
        if f == 0:
            break
        x = 1/f
    return out

t0 = time.time()
# ---- Step 1: enclosure from two convergents (k = 24 below, k = 23 above),
# found once by any means; the certificate only needs the two comparisons.
a_lo, b_lo = 1193652440098, 753110839881      # p_24 / q_24 < theta
a_hi, b_hi = 217976794617, 137528045312       # p_23 / q_23 > theta
ok1, enc1 = less_pow2_pow3(a_lo, b_lo)
ok2, enc2 = less_pow3_pow2(b_hi, a_hi)
lo, hi = Fr(a_lo, b_lo), Fr(a_hi, b_hi)
print('Step 1: 2^a < 3^b  for (a,b) =', (a_lo, b_lo), ':', ok1)
print('        3^b < 2^a  for (b,a) =', (b_hi, a_hi), ':', ok2)
print('        hi - lo =', float(hi - lo), ' (exactly 1/(b b\') =', hi - lo == Fr(1, b_lo*b_hi), ')')
assert ok1 and ok2

# ---- Step 2: ln 2 enclosure
N = 90
l2 = sum(Fr(1, n * 2**n) for n in range(1, N+1))
u2 = l2 + Fr(1, (N+1) * 2**N)
print('Step 2: ln 2 in [l2, u2], width', float(u2 - l2), ', l2 =', float(l2))
# short constants: 69/100 < ln 2 < 7/10, from seven terms of the same series
s7 = sum(Fr(1, n * 2**n) for n in range(1, 8))
assert Fr(69, 100) < s7 and s7 + Fr(1, 8 * 2**7) < Fr(7, 10)
print('        short form: sum_{n<=7} 1/(n 2^n) =', s7, '=', float(s7), '> 69/100; plus tail 1/1024 gives',
      float(s7 + Fr(1, 1024)), '< 7/10')

# ---- Step 3: partial quotients
cl, ch = cf_rational(lo), cf_rational(hi)
K = 0
while K < min(len(cl), len(ch)) and cl[K] == ch[K]:
    K += 1
A = cl[:K-1]           # drop the last common term to be safe
print('Step 3: certified partial quotients a_0..a_%d =' % (len(A)-1), A)
assert A[:23] == [1, 1, 1, 2, 2, 3, 1, 5, 2, 23, 2, 2, 1, 1, 55, 1, 4, 3, 1, 1, 15, 1, 9]

p = [0, 1]; q = [1, 0]
for ai in A:
    p.append(ai*p[-1] + p[-2]); q.append(ai*q[-1] + q[-2])
P_ = lambda k: p[k+2]; Q_ = lambda k: q[k+2]
# determinant identities (Lemma basic (ii))
assert all(P_(k)*Q_(k-1) - P_(k-1)*Q_(k) == (-1)**(k-1) for k in range(0, len(A)))

# ---- Step 4: rungs and exclusion checks
rungs = []
for k in range(0, len(A), 2):
    # when a_{k+1} is not certified we use only a_{k+1} >= 1, i.e. the rung j = 1
    for j in range(1, (A[k+1] if k+1 < len(A) else 1) + 1):
        rungs.append((Q_(k-1) + j*Q_(k), P_(k-1) + j*P_(k), k, j))
rungs.sort()
rstar = None; table = []; nchk = 0; maxbits = 0
for (u, P, k, j) in rungs:
    dlo = P - u*hi            # delta(u) > dlo
    dhi = P - u*lo            # delta(u) < dhi
    assert 0 < dlo < dhi < 1  # so ceil(u theta) = P exactly
    excl = dlo * 3 * X * l2 >= u      # rational inequality, i.e. integer cross multiplication
    num = (dlo * 3 * X * l2).numerator
    maxbits = max(maxbits, num.bit_length())
    nchk += 1
    table.append((u, P, k, j, dlo, dhi, excl))
    if not excl:
        feasible = dhi * 3 * X * u2 < u   # certifies delta(u) < u/(3 X ln 2)
        rstar = u
        print('Step 4: first rung not excluded: r_* =', u, ' s =', P, ' (k=%d, j=%d)' % (k, j),
              ' certified feasible:', feasible)
        print('        from the upper end of the enclosure: delta(r_*) < %.4e, margin delta(r_*)/(r_* eps) < %.4f'
              % (float(dhi), float(dhi*3*X*u2/u)))
        break
nr = len(table)
print('        rungs checked:', nr, ' all before r_* excluded:', all(t[6] for t in table[:-1]))
# margins: ratio delta_lo(u) / (u eps) using l2 (exact rationals, printed as floats)
marg = [(t[0], float(t[4]*3*X*l2/t[0])) for t in table[:-1]]
mm = min(marg, key=lambda z: z[1])
print('        smallest certified margin delta(u)/(u eps) =', mm[1], ' at u =', mm[0])
print('        largest integer in the exclusion checks: %d bits' % maxbits)

# ---- Step 4 again with the short constants 69/100 and 7/10 in place of l2, u2
ok_short = True
for (u, P, k, j, dlo, dhi, excl) in table[:-1]:
    if not dlo * 3 * X * Fr(69, 100) >= u: ok_short = False
u, P = table[-1][0], table[-1][1]
feas_short = (P - u*lo) * 3 * X * Fr(7, 10) < u
print('        with ln 2 replaced by 69/100 (exclusions) and 7/10 (feasibility):', ok_short, feas_short)
# integer form of one check, for display: rung q_21
u, P = 6586818670, 10439860591

num = (P*b_hi - u*a_hi) * 3 * 2**68 * 69     # (P - u a_hi/b_hi) 3 X (69/100) >= u  <=>  num >= 100 u b_hi
print('        integer check at q_21:', num, '>=', 100*u*b_hi, ':', num >= 100*u*b_hi,
      ' ratio', float(Fr(num, 100*u*b_hi)))

# ---- Step 5: consistency
dec = all(table[i+1][5] < table[i][4] for i in range(nr-1))
print('Step 5: delta enclosures strictly decreasing along the rungs:', dec)
SCAN = 2*10**6
recs = []; cur = None; undecided = 0
for r in range(1, SCAN+1):
    s_lo = -((-r*a_lo)//b_lo); s_hi = -((-r*a_hi)//b_hi)
    assert s_lo == s_hi            # ceil(r theta) certified
    dl, dh = s_lo - r*hi, s_lo - r*lo
    if cur is None or dh < cur[0]:
        recs.append(r); cur = (dl, dh)
    elif not dl > cur[1]:
        undecided += 1
print('        undecided comparisons in the scan:', undecided)
lad = [t[0] for t in rungs if t[0] <= SCAN]  # rungs below SCAN only need k <= 8
print('        exact scan of r <= %d: records = ladder rungs:' % SCAN, recs == lad, ' (%d records)' % len(recs))
print('total time %.2f s' % (time.time() - t0))
print('\nSelected rungs (u, s, k, j, certified bound for delta, margin): lower bounds for u < r_*,')
print('upper bounds (from the other end of the enclosure, with u2) at u = r_*')
for t in table:
    if t[0] >= 15601 and t[0] != rstar:
        print('  %16d %16d  k=%2d j=%2d  delta>%.6e  margin>%.4f' % (t[0], t[1], t[2], t[3], float(t[4]), float(t[4]*3*X*l2/t[0])))
    elif t[0] == rstar:
        print('  %16d %16d  k=%2d j=%2d  delta<%.6e  margin<%.4f' % (t[0], t[1], t[2], t[3], float(t[5]), float(t[5]*3*X*u2/t[0])))
print('enclosure of 3^b_lo: L, U mantissas of', enc1[0].bit_length(), 'bits, exponent', enc1[1], ', U - L =', enc1[2]-enc1[0] if enc1[1]==enc1[3] else 'n/a')
print('squarings per power:', len(bin(b_lo))-3, len(bin(b_hi))-3)
