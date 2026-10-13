# Joint spectral radius of the SMM step matrices on one eigendirection, a = alpha*lambda.
# D: descent step  e' = (1-a) e           state (e_t, e_{t-1})
# N: injection     e' = (1-a)((1+b) e - b e_prev)
# Lower bound: max over products of length <= K of rho(P)^(1/len); upper bound: min over k of
# max over products of length k of ||P||^(1/k) in an optimized ellipsoidal norm (from SDP free iteration).
import numpy as np, itertools, sys
b0 = np.log(1.5)/np.log(2)
def mats(a, b):
    D = np.array([[1-a, 0.0], [1.0, 0.0]])
    N = np.array([[(1-a)*(1+b), -(1-a)*b], [1.0, 0.0]])
    return D, N
def rho(P): return max(abs(np.linalg.eigvals(P)))
def jsr_bounds(a, b, K=14):
    D, N = mats(a, b); lo = 0.0; best_up = np.inf
    prods = {(): np.eye(2)}
    for k in range(1, K+1):
        new = {}
        for w, P in prods.items():
            for s, M in (('D', D), ('N', N)):
                Q = M @ P; new[w+(s,)] = Q
        prods = new
        lo = max(lo, max(rho(P)**(1.0/k) for P in prods.values()))
        up = max(np.linalg.norm(P, 2) for P in prods.values())**(1.0/k)
        best_up = min(best_up, up)
    return lo, best_up
if __name__ == '__main__':
    for b in (b0, 0.9):
        print(f"beta={b:.3f}")
        for a in (0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0):
            lo, up = jsr_bounds(a, b)
            D, N = mats(a, b)
            print(f"  a={a:5}: JSR in [{lo:.4f}, {up:.4f}]  rho(D)={rho(D):.4f} rho(N)={rho(N):.4f} "
                  f"per step rho(ND)^(1/2)={np.sqrt(rho(N@D)):.4f}  rho(NN)^(1/2)={rho(N):.4f}")

# Exact JSR via the rank one structure D = v u^T, u = (1,0), v = (1-a, 1):
# any cyclic product containing D has spectral radius prod_i |u^T N^{j_i} v|, so
# JSR = max( rho(N), sup_j |u^T N^j v|^{1/(j+1)} ).
def jsr_exact(a, b, J=4000):
    D, N = mats(a, b); u = np.array([1.0, 0.0]); v = np.array([1-a, 1.0])
    best = rho(N); arg = 'N'; w = v.copy()
    for j in range(J):
        c = abs(u @ w)
        if np.linalg.norm(w) < 1e-250: break     # avoid underflow artefacts
        val = c**(1.0/(j+1)) if c > 0 else 0.0
        if val > best + 1e-15: best, arg = val, j
        w = N @ w
    return best, arg
if __name__ == '__main__':
    print("exact JSR = max(rho(N), sup_j |u^T N^j v|^(1/(j+1)))")
    for b in (b0, 0.9, 0.99):
        worst_excess = -1; rows = []; worst = None
        for a in np.concatenate([np.geomspace(1e-4, 0.999, 400)]):
            J, arg = jsr_exact(a, b)
            D, N = mats(a, b)
            ref = max(rho(D), rho(N))
            if J - ref > worst_excess: worst_excess, worst = J - ref, (a, J, ref)
            if a in (1e-4,) or abs(np.log10(a) - round(np.log10(a))) < 1e-9: pass
            rows.append((a, J, arg, rho(D), rho(N)))
        print(f"beta={b:.3f}: max over a of JSR - max(rho(D),rho(N)) = {worst_excess:.3e}; max JSR = {max(r[1] for r in rows):.6f}")
        a_w, J_w, ref_w = worst     # cost of the largest excess per step and in iteration count
        print(f"   largest excess at a={a_w:.4f}: JSR={J_w:.5f} against {ref_w:.5f}; factor larger by "
              f"{100*(J_w/ref_w-1):.2f}%, iterations by {100*(np.log(ref_w)/np.log(J_w)-1):.2f}%")
        for a in (0.001, 0.01, 0.1, 0.3, 0.5, 0.8):
            J, arg = jsr_exact(a, b); D, N = mats(a, b)
            print(f"   a={a}: JSR={J:.6f} (attained by {'pure N' if arg=='N' else 'D N^%d' % arg})  rho(D)={rho(D):.6f} rho(N)={rho(N):.6f}")
