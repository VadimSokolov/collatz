import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
import numpy as np, sys
sys.path.insert(0, BOOK + '/scripts')
def syr(x):
    m = 3*x+1; k = (m & -m).bit_length()-1
    return m >> k
# exact restricted matrix
M = 27
units = [r for r in range(M) if r % 3]
order = 18; inv2 = pow(2, -1, M)
P = {}
for r in units:
    y = (3*r+1) % M
    for k in range(1, order+1):
        y = y*inv2 % M
        P[(y, r)] = P.get((y, r), 0) + 2.0**-k/(1-2.0**-order)
odd = [r for r in units if r % 2 == 1]
R = np.array([[P.get((s, r), 0) for r in odd] for s in odd])
R = R / R.sum(axis=0)
print("exact restricted 9 state:", np.round(np.sort(np.abs(np.linalg.eigvals(R)))[::-1][:5], 4))
# source procedure: until 400000 admissible transitions
idx = {r: i for i, r in enumerate(odd)}
C = np.zeros((9, 9)); drawn = 0; x = 1
while drawn < 400000:
    x += 2
    if x % 3 == 0: continue
    a, b = x % M, syr(x) % M
    if a in idx and b in idx:
        C[idx[b], idx[a]] += 1; drawn += 1
P9 = C/C.sum(axis=0)
print("source procedure: last x", x, np.round(np.sort(np.abs(np.linalg.eigvals(P9)))[::-1][:5], 4))
