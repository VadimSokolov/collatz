"""Table tab:ch9:moduli: model chains modulo 3^a (halving exponents independent Geometric(1/2)), a = 2,3,4,5.
(i) restricted chain on odd representatives r in [1,3^a) coprime to 3, as in
    the source computation: empirical over odd n <= 1e6 (pairs with both
    residues odd representatives) and exact model (restrict the exact model
    matrix and renormalise columns); report second largest eigenvalue modulus.
(ii) full chain on all units mod 3^a: exact model matrix P satisfies P^a of
    rank one; report max column TV between empirical and exact matrices."""
import numpy as np
def syr(x):
    m = 3*x+1; k = (m & -m).bit_length()-1
    return m >> k
NMAX = 10**6
pairs = [(x, syr(x)) for x in range(1, NMAX+1, 2) if x % 3]
def exact(M):
    units = [r for r in range(M) if r % 3]; id_ = {r: i for i, r in enumerate(units)}
    order = 2*M//3; inv2 = pow(2, -1, M); P = np.zeros((len(units),)*2)
    for r in units:
        y = (3*r+1) % M
        for k in range(1, order+1):
            y = y*inv2 % M; P[id_[y], id_[r]] += 2.0**-k/(1-2.0**-order)
    return units, id_, P
def lam2(P): return np.sort(np.abs(np.linalg.eigvals(P)))[::-1][1]
for a in [2, 3, 4, 5]:
    M = 3**a; units, id_, P = exact(M)
    odd = [r for r in units if r % 2]; io = {r: i for i, r in enumerate(odd)}
    R = P[np.ix_([id_[r] for r in odd], [id_[r] for r in odd])]; R = R/R.sum(axis=0)
    C = np.zeros((len(odd),)*2); Cf = np.zeros_like(P)
    for x, s in pairs:
        u, v = x % M, s % M
        Cf[id_[v], id_[u]] += 1
        if u % 2 and v % 2: C[io[v], io[u]] += 1
    Re = C/C.sum(axis=0); Pe = Cf/Cf.sum(axis=0)
    rank = np.linalg.matrix_rank(np.linalg.matrix_power(P, a), tol=1e-9)
    rank1 = np.linalg.matrix_rank(np.linalg.matrix_power(P, a-1), tol=1e-9)
    tv = 0.5*np.abs(Pe-P).sum(axis=0).max()
    print("M=%3d | restricted: states %2d emp |l2| %.4f exact |l2| %.4f | full: states %3d rank P^a %d rank P^(a-1) %d maxTV(emp,exact) %.4f min count %d"
          % (M, len(odd), lam2(Re), lam2(R), len(units), rank, rank1, tv, Cf.sum(axis=0).min()))
# Chernoff rate for the stopping time tail of the Bernoulli model on T
th = np.log(np.log(2)/np.log(1.5))/np.log(3)
rho = 0.5*(2**-th + 1.5**th)
print("Chernoff: theta_min %.4f  rho %.4f" % (th, rho))
