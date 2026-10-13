# Stopping times, Terras bijection, coefficient stopping time, Chernoff bound
import math
def T(n): return n//2 if n%2==0 else (3*n+1)//2
def parvec(n,k):
    v=[]
    for _ in range(k): v.append(n%2); n=T(n)
    return v
# Terras bijection check
for k in range(1,17):
    vs=set(tuple(parvec(r,k)) for r in range(2**k))
    assert len(vs)==2**k
print("bijection ok k<=16")
# affine formula check
import random
random.seed(1)
for _ in range(2000):
    n=random.randrange(1,10**9); k=random.randrange(1,40)
    x=n; m=0; rho=0
    for i in range(k):
        if x%2: rho+= 0  # placeholder
        x=T(x)
    v=parvec(n,k); idx=[i for i in range(k) if v[i]]; m=len(idx)
    rho=sum(3**(m-j)*2**idx[j-1] for j in range(1,m+1))
    assert 2**k*x==3**m*n+rho
print("affine ok")
def tau_vec(v):
    m=0
    for j,b in enumerate(v,1):
        m+=b
        if 3**m<2**j: return j
    return None
def sigma(n):
    x=n;k=0
    while True:
        x=T(x);k+=1
        if x<n: return k
# table of residues with tau = k
cum=0
rows=[]
for k in range(1,21):
    res=[r for r in range(2**k) if tau_vec(parvec(r,k))==k]
    cum+=len(res)/2**k
    rows.append((k,len(res),res if len(res)<=8 else None,1-cum))
    print(k,len(res),res if len(res)<=10 else '...', "P(tau>k)=%.6f"%(1-cum))
# sigma == tau for 2<=n<=10^6
N=10**6; bad=0
for n in range(2,N+1):
    s=sigma(n)
    # tau via parity vector up to s
    v=parvec(n,s)
    if tau_vec(v)!=s: bad+=1
print("sigma==tau failures for 2..1e6:",bad)
# stopping time distribution empirical for n<=1e6
from collections import Counter
c=Counter(sigma(n) for n in range(2,N+1))
print(sorted(c.items())[:20], "max sigma", max(c))
mx=max(range(2,N+1),key=sigma); print("argmax sigma <=1e6",mx,sigma(mx))
# Chernoff
th=math.log(2)/math.log(3)
D=th*math.log(2*th)+(1-th)*math.log(2*(1-th))
print("theta",th,"D",D)
for k in [10,20,40,100]:
    # exact binomial tail P(Bin(k,1/2) >= theta k)
    t=sum(math.comb(k,i) for i in range(math.ceil(th*k),k+1))/2**k
    print(k,"binom tail",t,"chernoff",math.exp(-k*D))
