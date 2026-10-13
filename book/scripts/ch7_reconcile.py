import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
# Observed ratios to GD in the noise experiment versus the predicted stationary ratios (see ch7_variance.py), and
# observed logistic speedups; the observations are read from ch7_experiments.json (script ch7_experiments.py).
import json
ex = json.load(open(BOOK + '/scripts/ch7_experiments.json'))
fin = ex['noise']['final']; gd = [m for m, s in fin['GD'][1:]]; smm = [m for m, s in fin['SMM'][1:]]
for k, v in fin.items():
    v = [m for m, s in v[1:]]
    r = [x / g for x, g in zip(v, gd)]
    print(k, 'ratio to GD', [round(x, 3) for x in r], 'range %.2f to %.2f' % (min(r), max(r)),
          'vs SMM', [round(x / s, 3) for x, s in zip(v, smm)])
for s in [0.05, 0.1, 0.2, 0.4]: print('GD predicted', s, 0.3054973 * s * s)
lg = {k: v[0] for k, v in ex['logistic']['counts'].items()}
for k, v in lg.items(): print(k, 'iterations %.1f' % v, 'speedup over GD %.3f' % (lg['GD'] / v))
for b in [0.5, 0.585, 2 / 3, 0.9]: print('underdamped threshold a >', 1 - 4 * b / (1 + b) ** 2)
