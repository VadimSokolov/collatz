"""Exact first nonexcluded r versus the crude bound min_j max(q_j, delta(q_j)/eps), for several X."""
from mpmath import mp, mpf, log, ceil, floor
mp.dps=80
th=log(3)/log(2)
a=[];x=th
for _ in range(40):
    ai=int(floor(x)); a.append(ai); x=1/(x-ai)
q=[0,1]
for ai in a[1:]: q.append(ai*q[-1]+q[-2])
c=set()
for k in range(1,len(q)-1):
    for j in range(0,a[k]+1): c.add(q[k-1]+j*q[k])
delta=lambda r: ceil(r*th)-r*th
recs=[];best=mpf(2)
for r in sorted(v for v in c if v>0):
    d=delta(r)
    if d<best: best=d; recs.append((r,d))
print("rungs (delta records) <= 72057431991:",sum(1 for r,_ in recs if r<=72057431991))
for e in [40,60,68,71,72,80]:
    X=mpf(2)**e; eps=1/(3*X*log(2))
    ex=next(r for r,d in recs if d<r*eps)
    cr=int(ceil(min(max(r,d/eps) for r,d in recs)))
    print(f"2^{e}: exact r = {ex}, s = {int(ceil(ex*th))}; crude = {cr}, crude s = {int(ceil(cr*th))}")
