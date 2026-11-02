"""Ladder of one sided best approximations s/r > theta = log2 3 and the
rational cycles of their Christoffel words (sharpness of Theorem 20.2)."""
import mpmath as mp, numpy as np, math
mp.mp.dps = 80
TH = mp.log(3, 2)
def cf(x, n):
    a = []
    for _ in range(n):
        q = int(mp.floor(x)); a.append(q); x = 1/(x-q)
    return a
A = cf(TH, 40)
# all intermediate fractions above theta, then keep records of delta(r)=ceil(r th)-r th
p0, q0, p1, q1 = 1, 0, A[0], 1
cands = set()
for k in range(1, len(A)):
    for j in range(0, A[k]+1):
        p, q = p0 + j*p1, q0 + j*q1
        if q > 0: cands.add((p, q))
    p0, q0, p1, q1 = p1, q1, A[k]*p1 + p0, A[k]*q1 + q0
rungs, best = [], mp.inf
for p, q in sorted(cands, key=lambda f: f[1]):
    if mp.mpf(p)/q <= TH: continue
    d = p - q*TH
    if p != int(mp.ceil(q*TH)): continue
    if d < best:
        best = d; rungs.append((q, p, d))
L2 = mp.log(2)
def F(r, d):  # crude bound r/(3 Lambda)
    return r/(3*d*L2)
if __name__ == "__main__":
    print("partial quotients:", A[:30])
    for r, s, d in rungs:
        if r > 2e13: break
        print(f"{r:>16d} {s:>16d} delta={mp.nstr(d,8):>14s} log2F={mp.nstr(mp.log(F(r,d),2),8)}")
