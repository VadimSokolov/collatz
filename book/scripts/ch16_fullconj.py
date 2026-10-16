"""P11: test the full conjecture  worst(N, beta) = max( 1/((4+2b)N+2), b^2/(4N-2+2(N-2)b) )  for even N,
and the thresholds beta_N = (sqrt(9N^2-2N+1) - N - 1)/(2N)."""
import sys, os, math, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ch5d_pep import pep, smm
rows = []
for N in (2, 4, 6, 8, 10):
    bN = (math.sqrt(9 * N * N - 2 * N + 1) - N - 1) / (2 * N)
    for beta in (0.3, 0.585, bN - 0.02, bN + 0.02, 0.95, 0.99):
        if beta >= 1: continue
        w = pep(smm(N, beta)); a = 1 / ((4 + 2 * beta) * N + 2); c = beta ** 2 / (4 * N - 2 + 2 * (N - 2) * beta)
        pred = max(a, c); rows.append((N, beta, w, pred, w / pred - 1))
        print(f"N={N:2d} beta={beta:.4f} (beta_N={bN:.5f}) pep={w:.8f} conj={pred:.8f} rel.err={w/pred-1:+.2e} [{'overshoot' if c > a else 'main'}]", flush=True)
json.dump(rows, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'ch16_fullconj.json'), 'w'))
