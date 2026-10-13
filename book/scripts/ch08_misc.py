"""Chapter 8 small checks: B3 range, crude vs exact bound components, rational cycle example."""
from fractions import Fraction
import mpmath as mp
def T(x): return x//2 if x%2==0 else (3*x+1)//2
mx=0;arg=0
for x in range(1,10001):
    y=x;j=0
    while y!=1: y=T(y);j+=1
    if j>mx: mx,arg=j,x
print("max T steps to reach 1 for 1<=x<=10^4:",mx,"at x =",arg)
# rational cycle with word (1,0,0)
def Tq(x):
    return x/2 if x.numerator%2==0 else (3*x+1)/2
x=Fraction(1,5); o=[x]
for _ in range(3): x=Tq(x); o.append(x)
print("orbit of 1/5:",o)
mp.mp.dps=60
TH=mp.log(3)/mp.log(2)
for e,q in [(40,190537),(68,6586818670)]:
    eps=1/(3*mp.mpf(2)**e*mp.log(2)); d=mp.ceil(q*TH)-q*TH
    print(e,q,"delta(q)/eps =",mp.nstr(d/eps,12))
print("3^53/4^42 =",mp.nstr(mp.mpf(3)**53/mp.mpf(4)**42,6),"3^10/4^8=",mp.nstr(mp.mpf(3)**10/4**8,4),"3^15/4^12=",mp.nstr(mp.mpf(3)**15/4**12,4))
print("1015/1458 =",mp.nstr(mp.mpf(1015)/1458,5), " 2^68 =", mp.nstr(mp.mpf(2)**68,6), " log2(2.1689e20)=",mp.nstr(mp.log(2.1689e20,2),5))
