from mpmath import mp, mpf, log, floor, ceil, power
from fractions import Fraction as F
import math, bisect, numpy as np
mp.prec=300
th=log(3)/log(2)
N=50000
d=[None]*(N+1)
for t in range(1,N+1):
    x=t*th; d[t]=ceil(x)-x
df=np.array([1.0]+[float(x) for x in d[1:]])
rungs=[];m=2
for r in range(1,N+1):
    if d[r]<m: m=d[r]; rungs.append(r)
# G'(r)=x0 of mechanical word K_t=floor(t th+delta(r)) : exact rational
def K(t,r): return int(floor(t*th+d[r])) if t>0 else 0
def x0(Ks,r,s):  # Ks[t], t=0..r-1
    c=sum(3**(r-1-j)*2**Ks[j] for j in range(r)); return F(c,2**s-3**r)
def cyc(Ks,r,s):
    x=x0(Ks,r,s);out=[x]
    for t in range(r-1):
        k=(Ks[t+1]-Ks[t]); x=(3*x+1)/2**k; out.append(x)
    return out
# exact check small r
ok_min=True;ok_cat=True;G={}
for r in range(1,400):
    s=int(ceil(r*th)); Ks=[K(t,r) for t in range(r)]
    c=cyc(Ks,r,s); G[r]=c[0]
    if min(c)!=c[0]: ok_min=False; 
    i=bisect.bisect_right(rungs,r)-1; ri=rungs[i]
    if ri!=r:
        u=r-ri; w=[K(t,u) for t in range(u)]+[int(ceil(u*th))+K(t,ri) for t in range(ri)]
        if w!=Ks: ok_cat=False
        assert G[r]<=max(G[ri],G[u])
print('x0 is min of mechanical cycle (r<400):',ok_min,' concatenation identity:',ok_cat)
best=-1;rec=[]
for r in range(1,400):
    if G[r]>best: best=G[r];rec.append(r)
print('exact records r<400',rec,[x for x in rungs if x<400])
print('ties at multiples', G[10]==G[5], G[15]==G[5], G[82]==G[41])
# book's R exact check at 10 vs r<10 in high precision
def Bk(r,dl):
    return 1+sum((power(2,d[t]-1) if d[t]>dl else 1) for t in range(1,r))
R=[None]+[Bk(r,d[r])/(3*log(2)*d[r]) for r in range(1,20)]
print('R(1..10)',[float(x) for x in R[1:11]])
print('R(10)>max R(<10):',R[10]>max(R[1:10]), float(R[10]/R[5]), float(1+(power(2,d[5])-1)*(Bk(5,d[5])-1)/(2*Bk(5,d[5]))))
# float version of G' to 50000: G'(r)=B'_r/(3(2^d-1)), B'=sum 2^{d-frac(t th+d)}
def Gp(r):
    dl=df[r]; dd=df[1:r]
    terms=np.where(dd>dl,2.0**(dd-1),2.0**(dd-1)*2.0)  # exceptional: 2^{delta(t)}
    return (1.0+terms.sum())/(3*(2.0**dl-1))
# check formula vs exact
print('float vs exact', max(abs(Gp(r)/float(G[r])-1) for r in range(1,400)))
best=-1;rec=[];vals={}
for r in range(1,20001):
    g=Gp(r);vals[r]=g
    if g>best*(1+1e-11): best=g;rec.append(r)
rs=[x for x in rungs if x<=20000]
print('float records to 20000 == rungs:',rec==rs)
print('book bound vs exact bound at rungs', [(r, float(Bk(r,d[r])/(3*log(2)*d[r])/G[r])) for r in rungs if r<400])
fails=[]
for r in range(1,400):
    s=int(ceil(r*th)); c=cyc([K(t,r) for t in range(r)],r,s)
    if min(c)!=c[0]: fails.append(r)
print('x0 not min at',fails[:30],len(fails)); print('fails at rungs:',[r for r in fails if r in rungs])
# general proof check: x_t>=x_0 iff ... check W_t<=0 for t>=1 at rungs
