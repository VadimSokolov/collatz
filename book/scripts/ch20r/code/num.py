import math
from scipy.optimize import brentq, minimize_scalar
def stats(m,d):
    mu=sum(math.log(x/d) for x in m)/d
    rho=lambda t: sum((x/d)**t for x in m)/d
    r=minimize_scalar(rho,bounds=(0,50),method='bounded')
    ts=None
    if mu<0:
        try: ts=brentq(lambda t:rho(t)-1,1e-6,200)
        except Exception as e: ts=float('inf')
    var=sum(math.log(x/d)**2 for x in m)/d-mu**2
    return mu,var,r.x,r.fun,-math.log(r.fun),-math.log2(r.fun),ts
for m,d in [((1,3),2),((1,5),2),((1,2,4),3),((1,2,7),3),((1,5,7),3),((1,3,5,7),4)]:
    print(m,d,stats(m,d), sum(m), d*d)
tc=math.log(2)/math.log(3); h=-tc*math.log2(tc)-(1-tc)*math.log2(1-tc)
print('tc',tc,'h',h,'1-h',1-h)
mu=0.5*math.log(.75); var=(math.log(3)**2)/4
print('mu',mu,'var',var,'2mu/var',2*abs(mu)/var,'1/mu',1/abs(mu))
eps=tc-.5; print('eps max',eps)
# Parrondo maps
A=lambda n: 3*n+1 if n%2==0 else n//4
B=lambda n: 3*n+1 if n%2==1 else n//4
def bounded(f,n,K=10000):
    s=set()
    for _ in range(K):
        if n in s: return True
        s.add(n); n=f(n)
    return False
print(all(bounded(A,n) for n in range(20000)), all(bounded(B,n) for n in range(20000)))
C=lambda n:B(A(n))
def lock(n):
    k=0
    while n%2: n=C(n);k+=1
    return k
print(max(lock(n) for n in range(200000)))
print([ (n,C(n)) for n in range(8)])
for f,nm in[(A,'A'),(B,'B'),(C,'C')]:
    n=27;o=[n]
    for _ in range(12): n=f(n);o.append(n)
    print(nm,o)
# cycles of A, B
def cyc(f,n):
    s=[]
    while n not in s: s.append(n); n=f(n)
    return tuple(sorted(s[s.index(n):]))
print(set(cyc(A,n) for n in range(5000)), set(cyc(B,n) for n in range(5000)))
# amusical
def g(n):
    if n%2==0: return 3*n//2
    if n%4==1: return (3*n+1)//4
    return (3*n-1)//4
for s in [1,2,4,44]:
    o=[s];n=g(s)
    while n!=s and len(o)<50: o.append(n);n=g(n)
    print(s,len(o),o)
n=8
for k in range(2000): n=g(n)
print(math.log10(n), 0.5*math.log(9/8)*2000/math.log(10))
# Terras bijection check d=3, m=(1,2,4)
def T(x,m=(1,2,4),d=3):
    i=x%d; r=(i*m[i])%d; return (m[i]*x-r)//d
k=5; seen=set()
for x in range(3**k):
    v=[];y=x
    for _ in range(k): v.append(y%3); y=T(y)
    seen.add(tuple(v))
print(len(seen),3**k)
# parrondo stationary
print(5/13*0.1+8/13*0.75)
# empirical stopping-time density tail for Collatz
import itertools
def frac(k):
    c=0
    for x in range(2**k):
        # residue class has stopping time >k iff 3^a/2^j >=1 for all j<=k
        y=x;a=0;ok=True
        for j in range(1,k+1):
            if y%2: a+=1; y=(3*y+1)//2
            else: y//=2
            if 3**a<2**j: ok=False;break
        c+=ok
    return c/2**k
for k in (10,16,20): 
    f=frac(k); print(k,f,-math.log2(f)/k)
