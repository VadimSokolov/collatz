"""chap:divergence: model chains on residues mod 3^a (halving exponents independent Geometric(1/2)).
(1) replicate the nine state computation (odd representatives mod 27);
(2) full chain on units mod 3^a: empirical (odd n <= 1e6) and exact model matrix;
(3) stationary law of the model chain (the Syracuse law, the law of Tao's Syracuse random
    variable) and fine scale comparison between scales."""
import numpy as np
from math import gcd

def syr(x):
    m = 3*x+1; k = (m & -m).bit_length()-1
    return m >> k, k

NMAX = 10**6
odds = range(3, NMAX+1, 2)
pairs = [(x, syr(x)[0]) for x in odds]

def eig_moduli(P):
    return np.sort(np.abs(np.linalg.eigvals(P)))[::-1]

# (1) nine state restricted chain, as in the source code (first 400000 admissible transitions, x not divisible by 3)
m = 27
states = [r for r in range(1, m, 2) if r % 3 != 0]
idx = {r: i for i, r in enumerate(states)}
C = np.zeros((9, 9)); drawn = 0
for x, s in pairs:
    if x % 3 == 0: continue
    a, b = x % m, s % m
    if a in idx and b in idx:
        C[idx[b], idx[a]] += 1; drawn += 1
        if drawn >= 400000: break
P9 = C / C.sum(axis=0)
print("(1) nine state chain, eigenvalue moduli:", np.round(eig_moduli(P9)[:5], 4), " transitions used", drawn)

def exact_matrix(M):
    units = [r for r in range(M) if r % 3 != 0]
    id_ = {r: i for i, r in enumerate(units)}
    order = 2*M//3  # 2 is a primitive root mod 3^a, order phi(M)
    P = np.zeros((len(units), len(units)))
    for r in units:
        a = (3*r+1) % M
        inv2 = pow(2, -1, M)
        y = a
        for k in range(1, order+1):
            y = (y*inv2) % M
            P[id_[y], id_[r]] += 2.0**(-k)/(1-2.0**(-order))
    return units, P

rows = []
for a in [2, 3, 4, 5]:
    M = 3**a
    units, Pex = exact_matrix(M)
    id_ = {r: i for i, r in enumerate(units)}
    Ce = np.zeros_like(Pex)
    for x, s in pairs:
        if x % 3 == 0: continue
        Ce[id_[s % M], id_[x % M]] += 1
    Pem = Ce / Ce.sum(axis=0)
    ee = eig_moduli(Pem); ex = eig_moduli(Pex)
    # exact mixing: P^a has identical columns
    Pa = np.linalg.matrix_power(Pex, a)
    spread = np.abs(Pa - Pa[:, [0]]).max()
    Pa1 = np.linalg.matrix_power(Pex, a-1)
    spread1 = np.abs(Pa1 - Pa1[:, [0]]).max()
    tvemp = 0.5*np.abs(Pem - Pex).sum(axis=0).max()
    print("(2) M=%d states %d | empirical |lambda_2| %.4f gap %.4f | exact |lambda_2| %.2e | max col diff P^a %.1e, P^(a-1) %.3f | max col TV(emp,exact) %.4f | min col count %d"
          % (M, len(units), ee[1], 1-ee[1], ex[1], spread, spread1, tvemp, Ce.sum(axis=0).min()))
    rows.append((M, len(units), ee[1], ex[1]))
    if a == 2:
        pi = Pa[:, 0]
        print("    Syracuse law mod 9:", {u: round(p*63, 4) for u, p in zip(units, pi)}, "(times 63)")

# (3) Syracuse laws pi_a and fine scale mixing TV(pi_a, uniform lift of pi_{a-1}), and TV(pi_a, uniform on units)
pis = {}
for a in range(1, 8):
    M = 3**a
    units, P = exact_matrix(M)
    v = np.zeros(len(units)); v[0] = 1
    for _ in range(a): v = P @ v
    pis[a] = dict(zip(units, v))
for a in range(2, 8):
    M = 3**a; Mp = 3**(a-1)
    tv_lift = 0.5*sum(abs(p - pis[a-1][u % Mp]/3) for u, p in pis[a].items())
    tv_unif = 0.5*sum(abs(p - 1/len(pis[a])) for p in pis[a].values())
    print("(3) a=%d  TV(pi_a, lift pi_{a-1}) %.4f   TV(pi_a, uniform units) %.4f" % (a, tv_lift, tv_unif))
tv12 = 0.5*sum(abs(p - pis[1][u % 3]/(3**6)) for u, p in pis[7].items())
print("(3) TV(pi_7, lift of pi_1) %.4f" % tv12)
print("pi_1:", pis[1])
