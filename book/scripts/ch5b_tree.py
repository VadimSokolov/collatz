"""Chapter 5b: doubling chain difference inequalities for the inverse tree.
Computes the certified growth exponent gamma_k for residue classes mod 3^k.
Operator: (Phi_g c)_r = min over the three lifts l of r to mod 3^(k+1) of
   sum_{j>=1, n 2^j = 2 or 8 mod 9} (3/2^(j+1))^g c[child(l,j)],
child(l,j) = (2^(j+1) l - 1)/3 mod 3^k.
gamma_k = sup{g : exists c>0 with Phi_g c >= c}.  Run: python3 -I ch5b_tree.py KMAX"""
import numpy as np, sys
from scipy.sparse import csr_matrix

def build(k):
    M = 3**k; M1 = 3**(k+1); P = 2*3**k            # period of 2^j mod 3^(k+1)
    classes = np.array([r for r in range(M) if r % 3])
    idx = -np.ones(M, dtype=int); idx[classes] = np.arange(len(classes))
    js = np.arange(1, P+1)
    pw1 = np.array([pow(2, int(j)+1, M1) for j in js], dtype=np.int64)
    pw9 = np.array([pow(2, int(j), 9) for j in js], dtype=np.int64)
    lifts = []
    for t in range(3):
        L = (classes + t*M).astype(np.int64)
        R, J = np.meshgrid(np.arange(len(classes)), js, indexing='ij')
        good = np.isin((pw9[None, :]*L[:, None]) % 9, (2, 8))
        child = (((pw1[None, :]*L[:, None] - 1) % M1)//3) % M
        rows = R[good]; jj = J[good]; cols = idx[child[good]]
        lifts.append((rows, cols, jj))
    return classes, lifts, P

def mats(lifts, n, P, g):
    tail = 1.0/(1.0 - 2.0**(-P*g))
    return [csr_matrix((tail*np.exp(g*(np.log(3.0)-(jj+1.0)*np.log(2.0))), (rows, cols)), shape=(n, n))
            for rows, cols, jj in lifts]

def certify(lifts, n, P, g, iters=600):
    W = mats(lifts, n, P, g)
    c = np.ones(n)
    for _ in range(iters):
        c2 = np.min([w @ c for w in W], axis=0); c = c2/c2.max()
    ratio = (np.min([w @ c for w in W], axis=0)/c).min()
    return ratio, c

def gamma_k(k, tol=2e-5):
    classes, lifts, P = build(k); n = len(classes)
    lo, hi = 0.05, 1.0
    while hi - lo > tol:
        mid = 0.5*(lo+hi)
        if certify(lifts, n, P, mid)[0] >= 1: lo = mid
        else: hi = mid
    return lo, classes, lifts, P

if __name__ == '__main__':
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    for k in range(1, kmax+1):          # k = 1 is the exercise; the table starts at k = 2
        g, classes, lifts, P = gamma_k(k)
        gc = np.floor(g*1000)/1000
        r, c = certify(lifts, len(classes), P, gc, 2000)
        print(f"k={k} classes={len(classes)} gamma_k~{g:.4f} certified at {gc:.3f}: "
              f"min ratio {r:.6f}, c range [{c.min():.4f},{c.max():.4f}]", flush=True)
    # exponent slightly above the tabulated one, used for the counting corollary
    for g in (0.47705, 0.4771):
        r, c = certify(lifts, len(classes), P, g, 2000)
        print(f"k={kmax} at g={g}: min ratio {r:.7f}", flush=True)
