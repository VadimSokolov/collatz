"""Numbers of the chapter The Price of a Quantifier from clock_sweep's output: quantiles, tail index, excursion fractions, correlation, stopping time records, total clock, tau(27), the Lagarias and Weiss constant, the one step factors of v^p*, and the two figures. Usage: python clock_stats.py sweep.bin figdir"""
import math, sys
import numpy as np
import matplotlib.pyplot as plt
d = np.fromfile(sys.argv[1]).reshape(-1, 4)
figdir = sys.argv[2]
N = len(d); n = np.arange(1, N + 1.0)
tau, M, sT, sO = d.T; sU = sT + sO
r = (tau / n**2)[1:]
print("quantiles", [round(float(np.quantile(r, q)), 3) for q in (.1, .25, .5, .75, .9, .99, .999)])
print("mean %.4g max %.4g at %d" % (r.mean(), r.max(), np.argmax(r) + 2))
for m in (10, 100, 1000):
    print("fraction t(n)>=%dn" % m, ((M[1:] / n[1:]) >= m).mean())
xs = np.sort(r)[::-1]; ys = np.arange(1, len(xs) + 1) / len(xs)
k = len(xs) // 10
print("tail index, upper decile", -np.polyfit(np.log(xs[:k]), np.log(ys[:k]), 1)[0])
for lo, hi in ((1e-1, 1e-2), (1e-2, 1e-3), (1e-3, 1e-4)):
    a, b = int(lo * len(xs)), int(hi * len(xs))
    print("local index", lo, hi, -np.polyfit(np.log(xs[b:a]), np.log(ys[b:a]), 1)[0])
print("correlation", np.corrcoef(np.log(r), 2 * np.log(M[1:] / n[1:]))[0, 1])
for lab, s in (("T0", sU), ("T", sT)):
    best = 0; out = []
    for i in range(1, N):
        v = s[i] / math.log(i + 1)
        if v > best:
            best = v; out.append((i + 1, int(s[i]), round(float(v), 2)))
    print("records", lab, out)
S = np.cumsum(tau)
for b in (10**3, 10**4, 10**5, 10**6):
    if b <= N:
        print("total clock N=%d %.4g, log_N %.3f" % (b, S[b - 1], math.log(S[b - 1]) / math.log(b)))
print("tau(27)", int(tau[26]), "t(27)", int(M[26]))
p = math.log(math.log(2) / math.log(1.5)) / math.log(3)
print("p*", p, "rho", .5 * (2**-p + 1.5**p))
for v in (1, 2, 3, 4, 5, 101, 10**4):
    print("one step factor at v=%d" % v, .5 * (.5**p + ((3 * v + 1) / (2 * v))**p))
f = lambda q: -q * math.log(q) - (1 - q) * math.log(1 - q) - q * math.log(3)
a, b = .55, .65
for _ in range(80):
    c = (a + b) / 2
    if f(a) * f(c) <= 0: b = c
    else: a = c
q = (a + b) / 2; g = 1 / (math.log(2) - q * math.log(3))
print("q", q, "gamma for T", g, "for T0", g * (1 + q))
plt.rcParams.update({"font.size": 9, "font.family": "serif"})
fig, ax = plt.subplots(figsize=(3.2, 2.6))
step = np.unique(np.logspace(0, math.log10(len(xs)), 2000).astype(int)) - 1
ax.loglog(xs[step], ys[step], color="black", lw=1, label=r"$\tau(n)/n^2$")
x0 = np.array([1e1, 1e9]); ax.loglog(x0, 0.8 * x0**-.5, color="gray", ls="--", lw=1, label="slope $-1/2$")
ax.set_xlabel(r"$x$"); ax.set_ylabel("fraction above $x$"); ax.legend(frameon=False)
fig.tight_layout(); fig.savefig(figdir + "/ch21_tail.pdf"); plt.close(fig)
fig, ax = plt.subplots(figsize=(3.2, 2.6))
ax.hist(np.log10(r), bins=80, color="0.45")
ax.set_xlabel(r"$\log_{10}\tau(n)/n^2$"); ax.set_ylabel("count")
fig.tight_layout(); fig.savefig(figdir + "/ch21_hist.pdf"); plt.close(fig)
