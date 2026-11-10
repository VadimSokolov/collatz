import numpy as np, math, io, contextlib
from exact_ratio import stat, beta
Vb=beta+2/(2+beta)
def lin(a):
    rho=(1-a)*(1-(1+beta)*a); ve=((1-a)**2*(1+beta)**2+1)/(1-rho**2); vo=(1-a)**2*ve+1
    return 0.5*(ve+vo)*a*(2-a)
for mu in (0.05,0.005):
    for s in (0.05,0.3,1.0,2.0,3.0,5.0):
        W=9*s*math.sqrt(1/mu)+40; dx=min(0.2,max(s/8,0.01))
        if s<=0.3: W=40*s+2; dx=s/40
        mg,ms,xs,n,bg,be=stat(mu,1.5,s,W=W,dx=dx)
        print(f'mu {mu} sigma {s} n {n}: GD {mg:.5g} SMM {ms:.5g} ratio {ms/mg:.4f}',flush=True)
    print('far field quadratic ratio',lin(mu),' V(beta)',Vb,'V^2',Vb**2)
