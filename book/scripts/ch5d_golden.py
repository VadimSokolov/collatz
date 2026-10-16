# Tight convex worst case of Switching MM (cyclic GD with steps h = 1+beta, 1) against the
# conjectured formula 1/((4+2 beta) N + 2) for even N, across beta, including beta beyond the
# golden threshold phi - 1 = 0.618.
import sys, os, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ch5d_pep import pep, smm
out = {}
for beta in (0.0, 0.3, 0.5, math.log(1.5)/math.log(2), 0.6, 0.618, 0.64, 0.7, 0.8, 0.9):
    rows = []
    for N in (2, 4, 8, 12, 16, 20, 24, 30):
        w = pep(smm(N, beta)); f = 1 / ((4 + 2 * beta) * N + 2)
        rows.append((N, w, f, w / f))
    out[f"{beta:.4f}"] = rows
    print(f"beta={beta:.4f} " + "  ".join(f"N={N}:{r:.4f}" for N, w, f, r in rows), flush=True)
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'ch5d_golden.json'), 'w'))
