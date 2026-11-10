import numpy as np, math, sys
from m4 import setup, Br
def dp(r,h,J,mode):
    m,g=setup(r); B=m.sum(); nb=int(B/h)+r+5
    a=np.zeros((J+1,nb)); a[0,0]=1.0
    for t in range(r):
        w=np.zeros_like(a)
        for j in range(J+1):
            val=m[t]*2.0**(-j)/h
            sh=int(math.floor(val)) if mode=='lo' else int(math.ceil(val))
            if sh==0: w[j]=a[j]
            else: w[j,sh:]=a[j,:-sh]
        if t==r-1: return w.sum(0),B
        suf=np.cumsum(w[::-1],0)[::-1]
        if g[t]==1: a=suf
        else:
            a=np.empty_like(w); a[0]=suf[0]; a[1:]=suf[:-1]
r=int(sys.argv[1]); h=1/float(sys.argv[2]); J=int(sys.argv[3])
lo,B=dp(r,h,J,'lo'); hi,_=dp(r,h,J,'hi')
print('r',r,'B',B,'total log10/r',math.log10(lo.sum())/r,math.log10(hi.sum())/r)
for c in [0.5,0.6,0.7,0.75,0.751,0.8,0.85,0.9,0.95,0.99]:
    k=int(math.ceil(c*B/h))
    nl=lo[k:].sum(); nh=hi[k:].sum()
    print(c, round(math.log10(nl)/r,4) if nl>0 else None, round(math.log10(nh)/r,4) if nh>0 else None)
