# Numbers for chapters 1 and 2
import math
from fractions import Fraction

def T0(n): return n//2 if n%2==0 else 3*n+1
def T(n): return n//2 if n%2==0 else (3*n+1)//2
def v2(n):
    c=0
    while n%2==0: n//=2; c+=1
    return c
def S(n):
    m=3*n+1; return m>>v2(m)

def orbit(n,f):
    o=[n]
    while o[-1]!=1: o.append(f(o[-1]))
    return o

o=orbit(27,T0)
print("27 T0 steps",len(o)-1,"max",max(o),"argmax step",o.index(max(o)),"odd steps",sum(x%2 for x in o[:-1]))
print("27 T steps",len(orbit(27,T))-1)
s=[27]
while s[-1]!=1: s.append(S(s[-1]))
print("27 Syracuse steps",len(s)-1)
print("orbit 27 first 30",o[:30])
for n in [3,6,7,9,15,97,871,6171,77031]:
    oo=orbit(n,T0); print(n,"steps",len(oo)-1,"max",max(oo))
print("orbit 7",orbit(7,T0))
print("orbit 7 T",orbit(7,T))
# records of total stopping time
def tst(n):
    c=0
    while n!=1: n=T0(n); c+=1
    return c
for N in [10**4,10**5,10**6]:
    best=max(range(1,N+1),key=tst); print("N",N,"argmax tst",best,tst(best))
# record holders (delay records) up to 1e5
rec=[];b=-1
for n in range(1,10**5+1):
    t=tst(n)
    if t>b: b=t; rec.append((n,t))
print("delay records",rec)
# highest peak below 1e4,1e5
for N in [10**4,10**5]:
    bm=0;bn=0
    for n in range(1,N+1):
        m=max(orbit(n,T0))
        if m>bm: bm=m;bn=n
    print("N",N,"max peak",bn,bm)
# geometric residues
for j in range(1,9):
    # v2(3n+1) >= j  <=> n = -3^{-1} mod 2^j
    c=(-pow(3,-1,2**j))%(2**j)
    # exactly j: mod 2^{j+1}
    e=[r for r in range(1,2**(j+1),2) if v2(3*r+1)==j]
    print("j",j,"ge class mod 2^j",c,"exact class mod 2^(j+1)",e)
# empirical distribution of k over odd n < 2^20
from collections import Counter
cnt=Counter(v2(3*n+1) for n in range(1,2**20,2))
print({k:cnt[k] for k in range(1,12)})
# 97 joins orbit of 27
o27=orbit(27,T0); o97=orbit(97,T0)
first=[x for x in o97 if x in set(o27)][0]
print("97 joins 27 orbit at",first,"after",o97.index(first),"steps; 27 reaches it after",o27.index(first))
print("27 check 3^41*27 < 2^70:",3**41*27<2**70, 3**41*27/2**70)
print("2^68 =",2**68)
# Syracuse orbit of 27, k values
s=[27];ks=[]
while s[-1]!=1:
    m=3*s[-1]+1;k=v2(m);ks.append(k);s.append(m>>k)
print("Syracuse orbit 27",s); print("k values",ks,"sum",sum(ks),"mean",sum(ks)/len(ks))
# negative cycles under T0 (python // and % work for negatives)
for n in [-1,-5,-17]:
    o=[n]; x=T0(n)
    while x!=n: o.append(x); x=T0(x)
    print("cycle from",n,"length",len(o),"min",min(o),"odd count",sum(1 for y in o if y%2))
print("peak 837799",max(orbit(837799,T0)))
# 5x+1 cycles among n<=100
def g(n): return n//2 if n%2==0 else 5*n+1
cyc=set()
for n in range(1,101):
    seen=[];x=n
    while x not in seen and x<10**12: seen.append(x); x=g(x)
    if x<10**12:
        i=seen.index(x); c=seen[i:]; cyc.add(min(c))
print("5x+1 cycle minima from n<=100:",sorted(cyc))
