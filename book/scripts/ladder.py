from mpmath import mp, mpf, log, ceil, floor
mp.dps=80
th=log(3)/log(2)
# continued fraction
a=[];x=th
for _ in range(40):
    ai=int(floor(x)); a.append(ai); x=1/(x-ai)
q=[0,1]  # q_{-1}=0,q_0=1
for ai in a[1:]: q.append(ai*q[-1]+q[-2])
cands=set()
for k in range(1,len(q)-1):
    for j in range(0,a[k]+1 if k<len(a) else 1):
        cands.add(q[k-1]+j*q[k])
cands=sorted(c for c in cands if c>0)
delta=lambda r: ceil(r*th)-r*th
recs=[];best=mpf(2)
for r in cands:
    d=delta(r)
    if d<best: best=d; recs.append((r,d))
print("a =",a[:25]); print("num delta records <2e13:", sum(1 for r,_ in recs if r<2e13))
# brute check records to 2e5
best=mpf(2);br=[]
mp.dps=40
for r in range(1,200001):
    d=delta(r)
    if d<best: best=d; br.append(r)
print("brute records match ladder up to 2e5:", br==[r for r,_ in recs if r<=200000])
mp.dps=80
for e,X in [(40,mpf(2)**40),(60,mpf(2)**60),(68,mpf(2)**68),(71,mpf(2)**71)]:
    eps=1/(3*X*log(2))
    exact=next(r for r,d in recs if d<r*eps)
    # note's crude bound: min_j max(q_j, delta(q_j)/eps) over record denominators
    crude=min(max(r, d/eps) for r,d in recs)
    print(f"2^{e}: exact first feasible r = {exact}, s=ceil(r th)={int(ceil(exact*th))} ; note crude bound = {int(ceil(crude))}")
r0=72057431991
print("bound at r0 = r/(3 Lambda) =", r0/(3*delta(r0)*log(2)), " log2 =", log(r0/(3*delta(r0)*log(2)))/log(2))
for e in [72,73,75,80]:
    X=mpf(2)**e; eps=1/(3*X*log(2))
    ex=next(r for r,d in recs if d<r*eps); print(f"2^{e}: r >= {ex}")
