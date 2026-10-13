"""Chapter 5b: the transient grid behind the transient bound conjecture.
tau_j = u^T N^j v is the first coordinate after j injection steps started from the state
v = (1-a, 1) left by one descent step: e_{t+1} = (1-a)((1+b) e_t - b e_{t-1}), e_0 = 1-a,
e_{-1} = 1.  On the grid of 200 values of b in [0, 0.999] and 300 values of a in
[1e-5, 0.9999] we iterate until a Lyapunov certificate shows that no later term can reach
1-a: with P solving P - N^T P N = I, the quantity V(z) = z^T P z does not increase along
z -> N z, and |u^T z|^2 <= (P^{-1})_{11} V(z), so (P^{-1})_{11} V(z_t) < (1-a)^2 bounds every
tau_j with j >= t.  Prints the number of certified pairs, the largest ratio
max_{j>=1} |tau_j|/(1-a) over the grid, and the number of violations.
Run: python3 -I ch5b_transient.py"""
import numpy as np

def lyapunov(a, b):
    """closed form solution P = [[p, q], [q, r]] of P - N^T P N = I for N = [[al, -ga], [1, 0]]"""
    al = (1 - a)*(1 + b); ga = (1 - a)*b
    p = 2/(1 - al**2 - ga**2 + 2*al**2*ga/(1 + ga))
    q = -ga*al*p/(1 + ga); r = 1 + ga**2*p
    return p, q, r

if __name__ == '__main__':
    B = np.linspace(0.0, 0.999, 200); A = np.geomspace(1e-5, 0.9999, 300)
    bb, aa = np.meshgrid(B, A, indexing='ij'); bb = bb.ravel(); aa = aa.ravel(); c = 1 - aa
    p, q, r = lyapunov(aa, bb); w = r/(p*r - q*q)            # w = (P^{-1})_{11}
    e = c.copy(); e_prev = np.ones_like(c); mx = np.zeros_like(c)
    active = np.ones(len(c), bool); steps = np.zeros(len(c), dtype=np.int64)
    t, CHUNK, TMAX = 0, 200, 20_000_000
    while active.any() and t < TMAX:
        idx = np.nonzero(active)[0]
        ci, bi, ei, epi, m = c[idx], bb[idx], e[idx], e_prev[idx], mx[idx]
        for _ in range(CHUNK):
            en = ci*(1 + bi)*ei - ci*bi*epi
            epi = ei; ei = en
            m = np.maximum(m, np.abs(en))                     # max over j >= 1 of |tau_j|
        t += CHUNK
        e[idx], e_prev[idx], mx[idx] = ei, epi, m
        V = p[idx]*ei*ei + 2*q[idx]*ei*epi + r[idx]*epi*epi
        ok = w[idx]*V < ci**2*(1 - 1e-9)
        active[idx[ok]] = False; steps[idx[ok]] = t
    ratio = mx/c; i = np.argmax(ratio)
    print(f"pairs {len(c)}, certified {int((~active).sum())}, largest number of steps needed {steps.max()}")
    print(f"max over the grid of max_(j>=1) |tau_j|/(1-a) = {ratio[i]:.7f} at beta={bb[i]:.4f}, a={aa[i]:.3g}")
    print(f"violations (ratio > 1): {int((ratio > 1).sum())}")
