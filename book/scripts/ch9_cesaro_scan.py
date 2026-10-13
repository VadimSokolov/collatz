"""Chapter 9: kbar at termination for odd 3 <= x0 <= 5001 (replicates O1 scan) and argmin,
and the largest rounding residual of the Cesaro identity eq:ch9:cesaro at termination."""
import math
import numpy as np
def syr(x):
    m = 3*x+1; k = (m & -m).bit_length()-1
    return m >> k, k
res = []; resid = 0.0
for x0 in range(3, 5002, 2):
    x, K, N, s = x0, 0, 0, 0.0
    while x != 1:
        s += math.log(3 + 1/x)
        x, k = syr(x); K += k; N += 1
    res.append((K/N, x0, K, N))
    resid = max(resid, abs(K*math.log(2) - math.log(x0) - s))
res.sort()
print(len(res), res[:3], np.mean([r[0] for r in res]))
print('largest residual of the Cesaro identity at termination: %.2e' % resid)
