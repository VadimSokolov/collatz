"""Ties Lambda_{tq} = t Lambda_q and how floating point breaks them (Chapter 8).

A tie F(tq) = F(q) is an exact equality between two irrational numbers, so at
any finite precision rounding decides it.  The script
  (1) prints the differences at the multiples 15, 918, 46803, 571611, 952685,
  (2) scans r <= 10^6 for strict records of F in double precision and at 150
      digits, with Lambda_r evaluated in the two forms
      s log 2 - r log 3 (s = floor(r theta) + 1) and (ceil(r theta) - r theta) log 2,
  (3) repeats the scan exactly: F(r) > F(r') if and only if
      r ceil(r' theta) > r' ceil(r theta), a comparison of integers once the
      ceilings are known (here from theta truncated to 10^-120).
"""
import math, mpmath as mp
mp.mp.dps = 150
TH = mp.log(3)/mp.log(2); L2m = mp.log(2); L3m = mp.log(3)
th = math.log(3)/math.log(2); L2 = math.log(2); L3 = math.log(3)
R = 10**6

def f1(r):  # s log 2 - r log 3 form, double precision
    s = math.floor(r*th) + 1; return r/(3*(s*L2 - r*L3))
def f2(r):  # (ceil(r theta) - r theta) log 2 form, double precision
    return r/(3*(math.ceil(r*th) - r*th)*L2)
def g1(r):  # s log 2 - r log 3 form, 150 digits
    s = int(mp.floor(r*TH)) + 1; return r/(3*(s*L2m - r*L3m))
def g2(r):  # (ceil(r theta) - r theta) log 2 form, 150 digits
    x = r*TH; return r/(3*(mp.ceil(x) - x)*L2m)
def fx(r):
    s = mp.ceil(r*TH); return r/(3*(s - r*TH)*mp.log(2))

for t, q in [(3, 5), (3, 306), (3, 15601), (3, 190537), (5, 190537)]:
    print(t*q, "double: form1 diff", f1(t*q) - f1(q), "form2 diff", f2(t*q) - f2(q),
          "| 150 digits: form1 diff", mp.nstr(g1(t*q) - g1(q), 3), "form2 diff", mp.nstr(g2(t*q) - g2(q), 3),
          "| delta(tq) - t delta(q) =", mp.nstr((mp.ceil(t*q*TH) - t*q*TH) - t*(mp.ceil(q*TH) - q*TH), 3))

A = [1, 1, 1, 2, 2, 3, 1, 5, 2, 23, 2, 2, 1, 1, 55, 1, 4, 3]
q = [0, 1]
for a in A[1:]: q.append(a*q[-1] + q[-2])
lad = set()
for k in range(0, len(q) - 2, 2):
    for j in range(1, A[k+1] + 1): lad.add(q[k] + j*q[k+1])
ladR = sorted(u for u in lad if u <= R)

def scan(F):
    run = 0; rec = []
    for r in range(1, R + 1):
        v = F(r)
        if v > run: run = v; rec.append(r)
    return rec

for name, F in [("double, form1", f1), ("double, form2", f2), ("150 digits, form1", g1), ("150 digits, form2", g2)]:
    rec = scan(F)
    print("%-18s strict records up to 10^6: %d; off the upper ladder: %s"
          % (name, len(rec), [r for r in rec if r not in lad]))

# exact scan: ceil(r theta) from theta truncated to 10^-120 (the truncation error,
# below 10^-114 for r <= 10^6, is far smaller than the least delta(r) here)
D = 10**120; T = int(mp.floor(TH*D))
best_r, best_c = 1, -((-T)//D)
rec = [1]
for r in range(2, R + 1):
    c = -((-r*T)//D)                  # ceil(r theta)
    if r*best_c > best_r*c:           # F(r) > F(best_r), exactly
        best_r, best_c = r, c; rec.append(r)
print("exact integer scan up to 10^6: %d strict records, equal to the upper ladder rungs: %s"
      % (len(rec), rec == ladR))
