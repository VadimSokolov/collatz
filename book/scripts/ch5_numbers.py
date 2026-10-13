"""Chapter 5 numbers: stride constants, orbit of 27 as overshoot and correction, identities."""
from math import log, sqrt
L = log(5/3); mu = abs(log(3/4)); k = L/mu
print("L", L, "mu", mu, "kappa", k, "rate", 1-k**-0.5, "opt coef", (sqrt(k)-1)/(sqrt(k)+1))
Lb = log(10/3); km = Lb/mu; print("mixed kappa", km, "rate", 1-km**-0.5)
Ls = log(2); ks = Ls/mu; print("kappa*", ks, "1-mu*/L*", 1-mu/Ls, "beta", log(1.5)/log(2), "1/(1-beta)", 1/(1-log(1.5)/log(2)))
b = log(1.5)/log(2)
print("((1+b)/(1-b))^2", ((1+b)/(1-b))**2, "(log3/log(4/3))^2", (log(3)/log(4/3))**2)
sig2 = 0.5*log(2)**2 + 0.5*log(1.5)**2 - (0.5*log(0.75))**2; print("sigma2", sig2, "C", sig2/(0.5*log(0.75))**2)
# orbit of 27 under standard map
n, odd, even, stridek, x = 27, 0, 0, [], 27
path = [n]
while n != 1:
    if n % 2: n = 3*n+1; odd += 1
    else: n //= 2; even += 1
    path.append(n)
print("27: steps", odd+even, "odd", odd, "even", even, "halvings per odd", even/odd, "max", max(path))
# strides of 27
y = 27; ks_ = []; ups = 0
while y != 1:
    m = 3*y+1; kk = 0
    while m % 2 == 0: m //= 2; kk += 1
    ks_.append(kk); ups += (kk == 1); y = m
print("strides", len(ks_), "ascending strides", ups, "kbar", sum(ks_)/len(ks_), "log2 3", log(3,2))
print("check cesaro", sum(ks_)*log(2), ">", len(ks_)*log(3)+log(27))
# longest run of ascending strides in orbit of 27
run = best = 0
for kk in ks_:
    run = run+1 if kk == 1 else 0; best = max(best, run)
print("longest ascending run", best, "k sequence", ks_)
# sensitivity of the stride supremum (rem:nesterovlimits): remove x = 3, then the limit of large x
def stride_inc(x):
    m = 3*x+1; kk = 0
    while m % 2 == 0: m //= 2; kk += 1
    return log(m) - log(x)
s5, a5 = max((stride_inc(x), x) for x in range(5, 10**6, 2))
print("sup over odd 5<=x<1e6", s5, "at", a5, "log(11/7)", log(11/7), "ratio", s5/mu, "rate", 1-(s5/mu)**-0.5)
print("large x limit log(3/2)/log(4/3)", log(1.5)/mu, "rate", 1-(log(1.5)/mu)**-0.5)
# largest one step ascent of T over odd x>=3 (tab:ch5:levels, row 2): log(5/3) at x = 3, as for strides
sT, aT = max((log((3*x+1)/2) - log(x), x) for x in range(3, 10**6, 2))
print("T step sup over odd 3<=x<1e6", sT, "at", aT, "log(5/3)", log(5/3), "kappa* at x=1", Ls/mu)
