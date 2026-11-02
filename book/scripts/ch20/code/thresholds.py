"""Verification range needed to exclude each ladder rung r, under four inequalities:
crude F(r) = r/(3 Lambda); run bound (r+s)/(9 Lambda); Hercher Thm 27, (3/4) F(r);
and Theorem 20.2, G(r) = B_r/(3 Lambda) with B_r = sum_{t<r} 2^{-{t theta}}."""
import numpy as np, mpmath as mp, math, re, os
from ladder import rungs, TH
mp.mp.dps = 50
def frac_parts(r):
    t = np.arange(r, dtype=np.float64)
    hi = float(TH - 1); lo = float((TH - 1) - mp.mpf(hi))
    prod = t*hi; ip = np.floor(prod); fr = (prod - ip) + t*lo
    return fr - np.floor(fr)
big = {}
fn = os.path.join(os.path.dirname(__file__), "../data/br_rungs.txt")
if os.path.exists(fn):
    for line in open(fn):
        m = re.search(r"R=(\d+)\s+B_lower=([\d.]+)\s+B_upper=([\d.]+)", line)
        if m: big[int(m.group(1))] = (float(m.group(2)) + float(m.group(3)))/2
rows = []
for r, s, d in rungs:
    if r > 1.4e11: break
    if r <= 2*10**7: B = float(np.exp2(-frac_parts(r)).sum())
    elif r in big: B = big[r]
    else: continue
    Lam = d*mp.log(2)
    F = r/(3*Lam)
    l2 = lambda v: float(mp.log(v, 2))
    rows.append((r, s, l2(F), l2((r+s)/(9*Lam)), l2(0.75*F), l2(B/(3*Lam)), B/r))
if __name__ == "__main__":
    print(f"{'r':>14} {'s':>14} {'log2 F':>9} {'run':>9} {'Hercher27':>9} {'Thm20.2':>9} {'B_r/r':>9}")
    for row in rows:
        print(f"{row[0]:>14d} {row[1]:>14d} {row[2]:9.4f} {row[3]:9.4f} {row[4]:9.4f} {row[5]:9.4f} {row[6]:9.6f}")
