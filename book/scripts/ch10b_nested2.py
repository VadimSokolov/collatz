import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
import numpy as np, math, json, time
tau = np.load(BOOK + '/scripts/data/ch3_tau_small.npy'); M = len(tau); B = 26
exact_max = int(tau[2:].max()); argmax = int(np.argmax(tau[2:]) + 2)
d = json.load(open(BOOK + '/scripts/data/ch10b_levels.json')); X = np.array(d['X']); D = math.log(2) + (math.log(2)/math.log(3))*math.log(math.log(2)/math.log(3)) + (1 - math.log(2)/math.log(3))*math.log(1 - math.log(2)/math.log(3))  # log 2 - H(1/log2 3)
k = np.arange(1, 401); lx = np.log(X)
print("local decay rates of X_k: [100,200] %.4f [200,300] %.4f [300,400] %.4f ; D = %.4f" % (
    -(lx[199]-lx[99])/100, -(lx[299]-lx[199])/100, -(lx[399]-lx[299])/100, D))
# fit log X_k + D k = a - p log k on k in [100,400]
A = np.vstack([np.ones(301), -np.log(k[99:])]).T; coef = np.linalg.lstsq(A, lx[99:] + D * k[99:], rcond=None)[0]
print("polynomial correction exponent p = %.3f (ballot type prediction 3/2)" % coef[1])
def nested(seed, nlive=500, steps=60, iters=5000, flips=(1, 2, 3)):
    rng = np.random.default_rng(seed)
    live = rng.integers(2, M, nlive); L = tau[live].astype(float); best = (0, 0); calls = 0; found = None
    for it in range(iters):
        i = int(np.argmin(L)); Lmin = L[i]
        j = i
        while j == i: j = int(rng.integers(nlive))
        x = int(live[j])
        for _ in range(steps):
            y = x
            for _ in range(int(rng.choice(flips))): y ^= (1 << int(rng.integers(0, B)))
            calls += 1
            if y >= 2 and tau[y] > Lmin: x = y
        live[i] = x; L[i] = tau[x]
        if L[i] > best[0]: best = (int(L[i]), x)
        if found is None and best[0] >= exact_max: found = (it, calls)
    return best, found, calls
def random_search(seed, calls):
    rng = np.random.default_rng(seed); xs = rng.integers(2, M, calls); return int(tau[xs].max())
for seed in range(3):
    best, found, calls = nested(seed)
    print(f"nested run {seed}: best tau {best[0]} at n={best[1]} with {calls} evaluations; exact max {exact_max} at {argmax}; found at {found}; random search same budget: {random_search(seed, calls)}")
