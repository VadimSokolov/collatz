import math, sys
import numpy as np

def T(v): return v//2 if v%2==0 else (3*v+1)//2
def c(v):
    if v%2==0:
        k=v//2; return 2*k*(k+1)
    k=(v-1)//2; return 5*k*k+9*k+4
def reset(n): return 2*n*(n+2)

def tau(n):
    s=0; v=n
    while v!=1:
        s+=c(v); v=T(v)
    return s

# 1. reconcile with the reported trace
total=0
for n in range(1,47):
    total+=tau(n)+reset(n)
print("formula time at completion of input 46:",total)
# walk input 47
t=total; v=47
while True:
    t+=c(v); nv=T(v)
    if nv==3644 or v==3644:
        print("value",v,"->",nv,"completes at",t)
    if v==3644: break
    v=nv

# 2. cost functional statistics
N=int(sys.argv[1]) if len(sys.argv)>1 else 10**6
ratio=np.empty(N); sig=np.empty(N,dtype=np.int64); mx=np.empty(N)
for n in range(1,N+1):
    s=0; v=n; steps=0; m=n
    while v!=1:
        s+=c(v); v=T(v); steps+=1
        if v>m: m=v
    ratio[n-1]=s/(n*n); sig[n-1]=steps; mx[n-1]=m
np.save("ratio.npy",ratio); np.save("sig.npy",sig); np.save("mx.npy",mx)
r=ratio[1:]
print("N",N)
for q in [0.1,0.25,0.5,0.75,0.9,0.99,0.999]:
    print("quantile",q,np.quantile(r,q))
print("mean",r.mean(),"max",r.max(),"argmax",np.argmax(ratio)+1)
# tail exponent: fit P(ratio>x) ~ x^-a for x in upper decile
xs=np.sort(r)[::-1]
k=len(xs)//10
y=np.log(np.arange(1,k+1)/len(xs)); x=np.log(xs[:k])
a=-np.polyfit(x,y,1)[0]
print("tail exponent estimate",a)
# 3. sigma records: U map stopping time = sigma_T + number of odd steps; compute sigma_U/ln n records
def sigmaU(n):
    s=0
    while n!=1:
        n = n//2 if n%2==0 else 3*n+1; s+=1
    return s
best=0
print("records of sigma_U(n)/ln n:")
for n in range(2,N+1):
    if n%2==0 and n<=N: pass
    s=sigmaU(n); rr=s/math.log(n)
    if rr>best:
        best=rr; print(n,s,round(rr,3))
