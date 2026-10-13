import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
# Figure ch2_stopdecay.pdf for Chapter 2 (shared book style): decay of the proportion
# 1 - F(k) of residues with coefficient stopping time tau > k, against the empirical
# fraction of n <= 10^6 with sigma(n) > k, the binomial tail of Step 3 and the Chernoff
# bound exp(-Dk) of thm:density. Replaces the ch2_stopdecay block of ch01_ch02_figures.py,
# without the band and the lower bound curve, which the text does not state.
import sys; sys.path.insert(0, BOOK + '/scripts'); from bookstyle import *
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

FIG = BOOK + "/figures/"
def T(n): return n // 2 if n % 2 == 0 else (3 * n + 1) // 2

K = 120
th = math.log(2) / math.log(3)
# exact 1 - F(k) by the recursion of prop:ch2:dp, with integer counts of alive words
st = {0: 1}; P = []
for j in range(1, K + 1):
    new = defaultdict(int)
    for m, c in st.items():
        for b in (0, 1):
            mm = m + b
            if 3 ** mm > 2 ** j: new[mm] += c
    st = new; P.append(sum(st.values()) / 2 ** j)
P = np.array(P)
Dk = th * math.log(2 * th) + (1 - th) * math.log(2 * (1 - th))
ks = np.arange(1, K + 1)
tail = np.array([sum(math.comb(k, i) for i in range(math.ceil(th * k), k + 1)) / 2 ** k
                 for k in map(int, ks)])
# stopping times sigma(n) of the accelerated map for 2 <= n <= 10^6
Nn = 10 ** 6
sig = np.zeros(Nn + 1, dtype=int)
for n in range(2, Nn + 1):
    m = n; c = 0
    while True:
        m = T(m); c += 1
        if m < n: break
    sig[n] = c
emp = np.array([(sig[2:] > k).mean() for k in ks])
sel = emp > 0
fig, ax = plt.subplots(figsize=(WIDTH, 2.9))
ax.semilogy(ks, np.exp(-ks * Dk), color=PAL['grey'], ls="--", lw=1.1, label="Chernoff bound $e^{-Dk}$")
ax.semilogy(ks, tail, color=PAL['gold'], lw=1.1, label=r"binomial tail $P(X\geq\theta k)$")
ax.semilogy(ks, P, color=PAL['teal'], lw=1.6, drawstyle="steps-post", label=r"exact $1-F(k)$")
ax.semilogy(ks[sel], emp[sel], color=PAL['coral'], lw=0, marker="o", ms=2.2,
            label=r"$n\leq 10^6$ with $\sigma(n)>k$")
ylo = 10 ** math.floor(math.log10(min(emp[sel].min(), P.min()) / 2))
ax.set_xlabel("$k$"); ax.set_ylabel("proportion"); ax.set_xlim(0, K); ax.set_ylim(ylo, 1.5)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), fontsize=7, ncol=2)
fig.savefig(FIG + "ch2_stopdecay.pdf"); plt.close(fig)
print("1-F(k) at k=20,40,100:", P[19], P[39], P[99])
print("ratio (1-F(k))exp(Dk) at k=20,100:", P[19] * math.exp(20 * Dk), P[99] * math.exp(100 * Dk))
print("D", Dk, "max sigma", sig.max(), int(np.argmax(sig)), "y lower limit", ylo)
print("fraction sigma=1, sigma=2:", (sig[2:] == 1).mean(), (sig[2:] == 2).mean())
