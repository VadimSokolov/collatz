import math
from collections import Counter
from scipy.optimize import brentq, minimize_scalar
A=lambda x: x//2 if x%2==0 else (3*x+1)//2
def mkB(m2):
    return lambda x: x//3 if x%3==0 else (x-1)//3 if x%3==1 else (m2*x-2)//3
def fate(f,n,K=20000,cap=10**80):
    seen=set()
    for _ in range(K):
        if n in seen: return 'cyc'
        if abs(n)>cap: return 'big'
        seen.add(n);n=f(n)
    return '?'
for m2 in (13,19):
    B=mkB(m2);C=lambda x:B(A(x))
    vals=[(0.5,1.5*m2/3),(1/6,1/6),(1/6,1/6),(1/6,m2/6)]
    g=lambda t:sum(p*v**(-t) for p,v in vals)
    ts=brentq(lambda t:g(t)-1,1e-4,20); r=minimize_scalar(g,bounds=(0,ts),method='bounded')
    print(m2,'theta*',round(ts,4),'rho',round(r.fun,5),'B',Counter(fate(mkB(m2),n) for n in range(1,20001)),'C',Counter(fate(C,n) for n in range(1,20001)))
    o=[27]
    for _ in range(10): o.append(C(o[-1]))
    print(o)
