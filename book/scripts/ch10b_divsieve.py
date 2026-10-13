# Divisibility sieve for cycles of the accelerated map: a cycle with m odd steps among k steps
# needs D = 2^k - 3^m to divide c_k(w) for its parity word w.  For a small prime p | D we compute,
# by dynamic programming over positions, the set of residues of c_k(w) mod p over all words with
# exactly m ones (and first letter 1, which loses no generality by rotation).  If 0 is not
# reachable, no cycle with these (m, k) exists.
import math
from sympy import primefactors, factorint
theta = math.log2(3)
def reachable_zero(m, k, p):
    # state: dict (ones, c mod p) -> True ; process positions t = 0..k-1 ; c_{t+1} = 3c + 2^t if odd
    S = {(1, 1 % p)}          # position 0 is odd: c_1 = 3*0 + 2^0 = 1
    for t in range(1, k):
        p2 = pow(2, t, p); new = set()
        for (o, c) in S:
            new.add((o, c))                                  # even step
            if o < m: new.add((o + 1, (3 * c + p2) % p))     # odd step
        S = {(o, c) for (o, c) in new if o + (k - 1 - t) >= m}   # prune hopeless
    return any(o == m and c == 0 for (o, c) in S), len(S)
excluded = []; tested = 0
for m in range(2, 41):
    k0 = math.ceil(m * theta)
    for k in range(k0, k0 + 4):
        D = 2 ** k - 3 ** m
        ps = [p for p in primefactors(D) if p <= 3000] if D < 10 ** 40 else []
        if not ps:
            # small factors only by trial division
            ps = [p for p in range(3, 3000, 2) if D % p == 0 and all(p % q for q in range(3, int(p ** 0.5) + 1, 2))]
        for p in ps:
            tested += 1
            ok, sz = reachable_zero(m, k, p)
            if not ok:
                excluded.append((m, k, p))
print("tested (m,k,p) triples:", tested, " excluded:", excluded[:20], len(excluded))
