"""Chapter 3 (O4): exact law of the first passage time of the Bernoulli walk
log x_n = u + sum Delta_i, Delta = log(3/2) or -log 2 with prob 1/2, below 0,
computed by dynamic programming over the number of up steps.  Returns mean,
variance, third and fourth cumulants for a given u, and block averages with
u distributed as log n for n uniform on [2^k, 2^{k+1}).  Output data/ch3_model.npz"""
import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
os.makedirs(BOOK + '/scripts/data', exist_ok=True)  # intermediate results that later scripts load
import numpy as np
L1, L2 = np.log(1.5), np.log(2.0)
m = 0.5*np.log(4/3); s2 = 0.25*np.log(3)**2

def passage_law(u, tol=1e-17):
    """P(tau = t) for t = 0,1,...; tau = first t with u + a L1 - (t-a) L2 <= 0."""
    alive = np.array([1.0]); probs = [0.0]; t = 0
    while True:
        t += 1
        new = np.zeros(t+1); new[:-1] += 0.5*alive; new[1:] += 0.5*alive
        a = np.arange(t+1); pos = u + a*L1 - (t-a)*L2
        dead = pos <= 0
        probs.append(new[dead].sum()); new[dead] = 0.0
        alive = new
        if t > 10 and alive.sum() < tol: break
        # trim leading zeros is not needed (vector small)
    return np.array(probs)

def cumulants(p):
    t = np.arange(len(p)); p = p/p.sum()
    mu = (t*p).sum(); d = t-mu
    k2 = (d**2*p).sum(); k3 = (d**3*p).sum(); k4 = (d**4*p).sum() - 3*k2**2
    return mu, k2, k3, k4

def block_law(k, npts=48):
    # n uniform on [2^k, 2^{k+1}) => log n has density e^u on [k L2, (k+1) L2]
    us = (k + (np.arange(npts)+0.5)/npts)*L2
    w = np.exp(us); w /= w.sum()
    laws = [passage_law(u) for u in us]
    T = max(len(l) for l in laws)
    P = np.zeros(T)
    for wi, l in zip(w, laws): P[:len(l)] += wi*l
    return P, (w*us).sum()

if __name__ == '__main__':
    import sys
    ks = [int(a) for a in sys.argv[1:]] or list(range(6, 30)) + [40, 50, 60, 100, 200, 500, 1000]
    rows = []
    for k in ks:
        P, Eu = block_law(k, 48 if k < 200 else 16)
        mu, k2, k3, k4 = cumulants(P)
        rows.append((k, Eu, mu, k2, k3, k4))
        print('k=%4d Eu %8.3f mean %9.3f (mean-u/m %.3f) var %10.2f (var-s2u/m3 %.2f) k3 %.0f (k3-3s4u/m5 %.0f) skew %.4f exkurt %.4f'
              % (k, Eu, mu, mu-Eu/m, k2, k2-s2*Eu/m**3, k3, k3-3*s2**2*Eu/m**5, k3/k2**1.5, k4/k2**2), flush=True)
    np.savez(BOOK + '/scripts/data/ch3_model.npz', rows=np.array(rows))
