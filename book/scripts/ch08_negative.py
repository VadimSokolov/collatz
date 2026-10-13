"""Cycles of the standard map T_0(n) = n/2 (n even), 3n+1 (n odd) on the negative integers."""
from fractions import Fraction
def T0(n): return n//2 if n%2==0 else 3*n+1
def cycle_from(n):
    seen={}
    i=0
    while n not in seen:
        seen[n]=i; n=T0(n); i+=1
    # extract cycle
    c=[n]; m=T0(n)
    while m!=n: c.append(m); m=T0(m)
    return tuple(c)
cycles={}
for n in range(-1,-10**6-1,-1):
    c=cycle_from(n); key=max(c)   # canonical: largest element (closest to 0)
    cycles.setdefault(key,(c,0)); cycles[key]=(cycles[key][0],cycles[key][1]+1)
for key,(c,cnt) in sorted(cycles.items(),reverse=True):
    # rotate to start at max
    i=c.index(key); c=c[i:]+c[:i]
    odd=[x for x in c if x%2]
    r=len(odd); s=len(c)-r
    prod=1
    for x in odd: prod*=Fraction(3*x+1,x)
    print("cycle",c); print("  length under T0",len(c),"odd elements",odd,"r=",r,"s=",s,"prod(3+1/x)=",prod,"2^s=",2**s,"starting values reaching it:",cnt)
print("positive 1..10^6 all reach cycle containing 1:", all(1 in cycle_from(n) for n in range(1,10**6+1)))
