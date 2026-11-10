import math, random, itertools
from fractions import Fraction as F
from scipy.optimize import minimize_scalar
A=lambda x: x//2 if x%2==0 else (3*x+1)//2
def mkB(m2):
    def B(x):
        i=x%3
        return x//3 if i==0 else (x-1)//3 if i==1 else (m2*x-2)//3
    return B
for m2 in (5,7,11,17,25):
    muB=math.log(m2/27)/3; muC=0.5*math.log(.75)+0.5*math.log(m2/3)+0.5*muB
    print(m2,'muB',round(muB,5),'muC',round(muC,6))
m2=7;B=mkB(m2);C=lambda x:B(A(x))
# exact symbol law over residues mod 2^k 3^k
k=4;N=2**k*3**k;cnt={}
for x in range(N):
    y=x+N*10**6  # any rep
    s=[]
    for _ in range(k):
        d=y%2;z=A(y);i=z%3;y=B(z);s.append((d,i))
    cnt[tuple(s)]=cnt.get(tuple(s),0)+1
p={(1,2):F(1,2),(0,0):F(1,6),(0,1):F(1,6),(0,2):F(1,6)}
ok=all(F(c,N)==math.prod(p.get(t,0) for t in s) for s,c in cnt.items())
print('iid law exact on residues mod 6^%d:'%k,ok,len(cnt),4**k)
# check symbols depend only on residue
import random
random.seed(1)
bad=0
for _ in range(2000):
    x=random.randrange(N);y1=x+N*random.randrange(10**9);y2=x+N*random.randrange(10**9)
    for _ in range(k):
        z1=A(y1);z2=A(y2)
        if (y1%2,z1%3)!=(y2%2,z2%3): bad+=1
        y1=B(z1);y2=B(z2)
print('bad',bad)
# simulation growth
tot=0;R=20000;n=0
for t in range(200):
    x=random.getrandbits(4000);l0=math.log(x)
    for _ in range(2000): x=C(x)
    tot+=(math.log(x)-l0)/2000
print('sim drift',tot/200)
# Chernoff
vals=[(0.5,3.5),(1/3,1/6),(1/6,7/6)]
r=minimize_scalar(lambda t:sum(p*v**(-t) for p,v in vals),bounds=(0,5),method='bounded');print('rho',r.x,r.fun)
var=sum(p*math.log(v)**2 for p,v in vals)-sum(p*math.log(v) for p,v in vals)**2;print('var',var)
# small integer behaviour
def fate(f,n,K=5000,cap=10**60):
    seen=set()
    for _ in range(K):
        if n in seen: return 'cyc'
        if abs(n)>cap: return 'big'
        seen.add(n);n=f(n)
    return '?'
from collections import Counter
print('B',Counter(fate(B,n) for n in range(1,20001)))
print('C',Counter(fate(C,n) for n in range(1,20001)))
def cyc(f,n):
    s=[]
    while n not in s:
        s.append(n);n=f(n)
        if abs(n)>10**60:return None
    return tuple(sorted(s[s.index(n):]))
print('Bcycles',set(cyc(B,n) for n in range(1,3000)))
cs=set(cyc(C,n) for n in range(1,3000));print('Ccycles',[c if c is None or len(c)<8 else (len(c),min(c)) for c in cs])
print([n for n in range(1,200) if fate(C,n)=='cyc'][:40])
