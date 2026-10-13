"""Numbers quoted in Chapters 4 and 5 (computed for the book)."""
import numpy as np
from math import log, sqrt
L2, L3 = log(2), log(3)

def syr(x):
    m = 3*x+1; k = 0
    while m % 2 == 0: m //= 2; k += 1
    return m, k

def v2(n):
    k=0
    while n%2==0: n//=2; k+=1
    return k

# initial ascent run = v2(x+1)-1, check
for x in range(3, 200001, 2):
    r, y = 0, x
    while True:
        y2, k = syr(y)
        if k != 1: break
        r += 1; y = y2
    assert r == v2(x+1)-1, x
print("initial ascent run = v2(x+1)-1 verified for odd x < 2e5")

# two stride: descent iff k1+k2>=4 ; density
def twostride(lo, hi):
    """Scan odd x with lo <= x <= hi."""
    cnt = 0; tot = 0; mx = -9; argmx = None; s = 0.0
    for x in range(lo, hi + 1, 2):
        y, k1 = syr(x); z, k2 = syr(y)
        d = log(z) - log(x)
        assert (d < 0) == (k1+k2 >= 4), x
        tot += 1; cnt += d < 0; s += d
        if d > mx: mx, argmx = d, x
    print("two stride, odd %d <= x <= %d: %d of %d descend, frac %.7f, mean %.6f, max %.4f at x=%d"
          % (lo, hi, cnt, tot, cnt/tot, s/tot, mx, argmx))
twostride(3, 20001)          # the 10000 starting values quoted in sec:ch4:twostride
N = 2*10**6
twostride(3, N - 1)          # all odd x with 3 <= x < 2*10^6

# log squared model
EDelta = log(3/4); VarDelta = 2*L2**2; ED2 = VarDelta + EDelta**2
V0 = ED2/(2*abs(EDelta))
print("E Delta", EDelta, "Var", VarDelta, "E D^2", ED2, "threshold logx", V0, "x", np.exp(V0))
# J strides
for J in (1,2):
    m = J*EDelta; v = J*VarDelta; e2 = v + m*m
    print("J",J,"threshold logx", e2/(2*abs(m)), "x", np.exp(e2/(2*abs(m))))

# MGF of stride increment
M = lambda p: 3**p/(2**(1+p)-1)
for p in (0.25, 0.5, 0.75, 1.0):
    print("p", p, "M", M(p), "(M-1)/p", (M(p)-1)/p)
from scipy.optimize import brentq
print("Wald root other than 0:", brentq(lambda p: M(p)-1, 0.5, 1.5))

# Ch 5 constants
kap = log(5/3)/abs(log(3/4)); print("kappa", kap, "rate", 1-kap**-0.5)
kap2 = log(10/3)/abs(log(3/4)); print("kappa mixed", kap2, "rate", 1-kap2**-0.5)
ks = L2/abs(log(3/4)); print("kappa*", ks, "rate", 1-ks**-0.5, "nesterov coef", (sqrt(ks)-1)/(sqrt(ks)+1), "1-1/k*", 1-1/ks)
beta = log(1.5)/L2; print("beta", beta, "log2 3 -1", L3/L2-1)
for b in (beta, 0.9, 0.5, 2/3):
    print("crossing beta", b, "kappa", ((1+b)/(1-b))**2)
for k in (kap, ks, kap2):
    print("opt coef at", k, (sqrt(k)-1)/(sqrt(k)+1))
print("KL U1||U2", L2, "KL U2||U3", 2*0.5*log(1.5))
