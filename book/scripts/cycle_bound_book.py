# Exact lower bound on the number r of odd elements of a nontrivial Syracuse
# cycle, given that every n <= X reaches 1.  By the mediant lemma only rungs of
# the one sided ladder of theta = log_2 3 need to be checked.
from mpmath import mp, mpf, log, ceil, floor
mp.dps = 80
theta = log(3) / log(2)

def partial_quotients(x, n=40):
    a = []
    for _ in range(n):
        ai = int(floor(x)); a.append(ai); x = 1 / (x - ai)
    return a

a = partial_quotients(theta)
q = [0, 1]                      # q_{-1} = 0, q_0 = 1
for ai in a[1:]:
    q.append(ai * q[-1] + q[-2])

# candidate rungs q_{k-1} + j q_k, then keep the records of delta
cands = sorted({q[k-1] + j * q[k] for k in range(1, len(q) - 1)
                for j in range(a[k] + 1)} - {0})
delta = lambda r: ceil(r * theta) - r * theta
rungs, best = [], mpf(2)
for r in cands:
    d = delta(r)
    if d < best:
        best = d; rungs.append((r, d))

def first_not_excluded(X):
    eps = 1 / (3 * X * log(2))
    return next(r for r, d in rungs if d < r * eps)

for e in (40, 60, 68, 71, 72):
    r = first_not_excluded(mpf(2) ** e)
    print(f"X = 2^{e}:  r >= {r},  s >= {int(ceil(r * theta))}")
