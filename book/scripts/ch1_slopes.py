# Chapter 1: least squares slopes, through the origin, of the total stopping
# time against log n, for the standard map T0 and the accelerated map T, over
# 2 <= n <= 10^4 (the range of the fit quoted in sec:ch1:orbits) and over
# 2 <= n <= 10^5 (the range of fig:ch1:totalstop). Heuristic slopes:
# 3/log(4/3) for T0 and 2/log(4/3) = 1/|log(sqrt(3)/2)| for T.
import math
import numpy as np

N = 10**5
def total_stopping_times(accelerated):
    st = np.zeros(N + 1, dtype=np.int64)
    for n in range(2, N + 1):
        x, c = n, 0
        while x >= n:            # run until the orbit drops below n
            if x % 2 == 0:
                x //= 2
            else:
                x = (3 * x + 1) // 2 if accelerated else 3 * x + 1
            c += 1
        st[n] = c + st[x]        # x < n is already known; st[1] = 0
    return st

for name, acc in (("standard T0", False), ("accelerated T", True)):
    st = total_stopping_times(acc)
    for top in (10**4, 10**5):
        ns = np.arange(2, top + 1)
        L, y = np.log(ns), st[2:top + 1].astype(float)
        print(f"{name:14s} 2<=n<={top:6d}  slope through origin {(L @ y) / (L @ L):.4f}")
    print(f"{name:14s} checks: sigma(27)={st[27]}, sigma(97)={st[97]}, max below 1e5 at n={int(np.argmax(st))} ({st.max()} steps)")
print(f"heuristic slopes: T0 {3 / math.log(4 / 3):.4f}, T {2 / math.log(4 / 3):.4f}")
