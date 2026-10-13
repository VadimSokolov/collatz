"""Chapter 9: Cesaro identity checks on long orbits, odd frequency thresholds,
entropy of the atypical parity frequency."""
import numpy as np
from math import log, log2
def syr(x):
    m = 3*x+1; k = (m & -m).bit_length()-1
    return m >> k, k
rng = np.random.default_rng(11)   # same seed and draws as ch9_figures.py
for bits in [1000, 2000, 4000]:
    x0 = int.from_bytes(rng.bytes(bits//8), 'big') | (1 << (bits-1)) | 1
    x = x0; N = 0; K = 0; corr = 0.0
    while x != 1:
        corr += log2(1 + 1/(3*x)); x, k = syr(x); N += 1; K += k
    lx = log2(x0)
    print("bits %d: N %d kbar %.5f  log2 3 + log2 x0/N = %.5f  correction/N = %.5f  identity residual %.1e"
          % (bits, N, K/N, log2(3) + lx/N, corr/N, K/N - (log2(3) + (lx + corr)/N)))
q = 1/log2(3); H = -q*log2(q) - (1-q)*log2(1-q)
print("odd freq threshold T: %.4f ; T0: %.4f ; H2(q) %.4f ; 1-H %.4f" % (q, 1/(1+log2(3)), H, 1-H))
