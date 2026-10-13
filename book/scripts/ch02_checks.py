# Checks for chapter 2 statements
from collections import Counter, defaultdict
def v2(n):
    c=0
    while n%2==0: n//=2; c+=1
    return c
cnt=Counter(v2(3*n+1) for n in range(1,2**20,2))
print(all(cnt[j]==2**(19-j) for j in range(1,20)), sum(c for j,c in cnt.items() if j>=20), sorted((j,c) for j,c in cnt.items() if j>=18))
def T0(n): return n//2 if n%2==0 else 3*n+1
for k in range(1,13):
    M=2**(k+2)
    groups=defaultdict(set)
    for n in range(M):
        w=[];x=n
        for _ in range(k): w.append(x%2); x=T0(x)
        groups[tuple(w)].add(n)
    ok=True
    for w,S in groups.items():
        if any(w[i]==1 and w[i+1]==1 for i in range(k-1)): ok=False
        z=w[:-1].count(0); mod=2**(z+1)
        res={n%mod for n in S}
        if len(res)!=1 or len(S)!=M//mod: ok=False
    fib=[0,1]
    while len(fib)<k+3: fib.append(fib[-1]+fib[-2])
    print(k,len(groups)==fib[k+2],ok)
