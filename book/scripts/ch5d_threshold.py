import sys, os, math, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ch5d_pep import pep, smm
def ratio(beta, N): return pep(smm(N, beta)) * ((4 + 2 * beta) * N + 2)
res = {}
for N in (2, 4, 6):
    lo, hi = 0.6, 0.999
    if ratio(hi, N) < 1 + 1e-5: res[N] = None; print(N, 'no break below 0.999'); continue
    for _ in range(22):
        mid = (lo + hi) / 2
        if ratio(mid, N) > 1 + 1e-5: hi = mid
        else: lo = mid
    res[N] = (lo + hi) / 2; print(N, 'threshold beta =', res[N], flush=True)
# closed form guess for N=2: compare the competing expression (1-h)^2-type tail
json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'ch5d_threshold.json'), 'w'))
