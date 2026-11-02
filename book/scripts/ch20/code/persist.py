"""Proposition 20.13: P(walk with steps log(3/2), -log 2 stays positive for k steps)
= P(u_j > j / theta for 1 <= j <= k), u_j = number of up steps; exact DP in floats."""
import numpy as np, math
TH = math.log2(3); D = -math.log(min(0.5*1.5**l + 0.5*2**-l for l in np.linspace(0.4,0.6,200001)))
K = 6000
p = np.zeros(K+2); p[0] = 1.0
out = {}
for j in range(1, K+1):
    q = np.zeros(K+2)
    q[1:] += 0.5*p[:-1]; q += 0.5*p          # up or down
    umin = math.floor(j/TH) + 1               # need u_j > j/theta
    q[:umin] = 0.0
    p = q
    if j in (100, 200, 400, 800, 1600, 3200, 6000):
        tot = p.sum(); out[j] = tot*math.exp(D*j)*j**1.5
print("D =", D)
for j, v in out.items(): print(j, round(v, 3))
