"""Chapter 5: capital dependent Parrondo games (Harmer and Abbott).
Exact long run drift (expected gain per game) for A, B, random mixture and periodic patterns."""
import numpy as np

def trans(game, eps):
    """Transition matrix on capital mod 3 and win probability vector."""
    if game == "A":
        p = np.full(3, 0.5 - eps)
    else:
        p = np.array([0.1 - eps, 0.75 - eps, 0.75 - eps])
    P = np.zeros((3, 3))
    for s in range(3):
        P[s, (s + 1) % 3] += p[s]
        P[s, (s - 1) % 3] += 1 - p[s]
    return P, p

def stationary(P):
    w, v = np.linalg.eig(P.T)
    pi = np.real(v[:, np.argmin(abs(w - 1))])
    return pi / pi.sum()

def drift_single(game, eps):
    P, p = trans(game, eps)
    pi = stationary(P)
    return pi, float(pi @ (2 * p - 1))

def drift_mix(eps, g=0.5):
    PA, pA = trans("A", eps); PB, pB = trans("B", eps)
    P = g * PA + (1 - g) * PB; p = g * pA + (1 - g) * pB
    pi = stationary(P)
    return pi, float(pi @ (2 * p - 1))

def drift_periodic(pattern, eps):
    mats = [trans(c, eps) for c in pattern]
    M = np.eye(3)
    for P, _ in mats:
        M = M @ P
    pi = stationary(M)  # distribution at start of each period
    tot, d = 0.0, pi.copy()
    for P, p in mats:
        tot += d @ (2 * p - 1)
        d = d @ P
    return float(tot / len(pattern))

for eps in (0.0, 0.005):
    print("eps", eps)
    print("  A drift", drift_single("A", eps)[1])
    pi, d = drift_single("B", eps); print("  B pi", np.round(pi, 4), "drift", round(d, 5))
    pi, d = drift_mix(eps); print("  random 1/2 mix pi", np.round(pi, 4), "drift", round(d, 5))
    for pat in ("AB", "AABB", "ABB", "AAB"):
        print("  periodic", pat, round(drift_periodic(pat, eps), 5))

# simulation check of AABB at eps 0.005
rng = np.random.default_rng(1)
cap, n = 0, 2_000_000
pat = "AABB"
for t in range(n):
    g = pat[t % 4]
    p = 0.5 - 0.005 if g == "A" else (0.1 - 0.005 if cap % 3 == 0 else 0.75 - 0.005)
    cap += 1 if rng.random() < p else -1
print("sim AABB drift", cap / n)

# exact rational checks at eps = 0
from fractions import Fraction as F
def exact_drift(p0, p1):
    # 3 state cycle chain: up with prob p_s, down with 1-p_s; p1 = p2
    import sympy as sp
    a, b, c = sp.symbols("a b c")
    q0, q1 = 1 - p0, 1 - p1
    eqs = [sp.Eq(a, c * p1 + b * q1), sp.Eq(b, a * p0 + c * q1), sp.Eq(a + b + c, 1)]
    s = sp.solve(eqs, [a, b, c])
    return s[a], s[a] * (2 * p0 - 1) + (1 - s[a]) * (2 * p1 - 1)
print("exact B", exact_drift(F(1, 10), F(3, 4)))
print("exact mix", exact_drift(F(3, 10), F(5, 8)))
print("B drift under uniform residues", F(1, 3) * (2 * F(1, 10) - 1) + F(2, 3) * (2 * F(3, 4) - 1))
