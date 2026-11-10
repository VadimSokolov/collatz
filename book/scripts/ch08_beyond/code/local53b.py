import numpy as np, math, sys
beta=math.log(1.5)/math.log(2); Vb=beta+2/(2+beta)
def run(g,xs,alpha,s,nch=4000,T=40000,burn=15000,seed=0):
    rng=np.random.default_rng(seed)
    x=np.full(nch,xs); acc=0;n=0
    for t in range(T):
        x=x-alpha*(g(x)+s*rng.standard_normal(nch))
        if t>=burn: acc+=((x-xs)**2).mean(); n+=1
    mgd=acc/n
    x=np.full(nch,xs); xp=x.copy(); acc=0; n=0
    for t in range(T):
        if t%2==1:
            y=x+beta*(x-xp); xn=y-alpha*(g(y)+s*rng.standard_normal(nch))
        else: xn=x-alpha*(g(x)+s*rng.standard_normal(nch))
        xp,x=x,xn
        if t>=burn: acc+=((x-xs)**2).mean(); n+=1
    return mgd,acc/n
def lin(a,alpha,s):
    rho=(1-a)*(1-(1+beta)*a); ve=alpha**2*s**2*((1-a)**2*(1+beta)**2+1)/(1-rho**2); vo=(1-a)**2*ve+alpha**2*s**2
    return alpha**2*s**2/(a*(2-a)), 0.5*(ve+vo)
L=1.0
for mu in [0.05,0.2]:
    g=lambda x: L*x-(L-mu)*np.tanh(x)      # f=L x^2/2-(L-mu)log cosh x : curvature mu at 0, L far away
    alpha=1/L; a=alpha*mu; vg,vs=lin(a,alpha,1.0)
    print('mu',mu,'linear ratio',vs/vg,'V(beta)',Vb)
    for s in [0.02,0.1,0.3,1.0,3.0,10.0]:
        mg,ms=run(g,0.0,alpha,s)
        print(f'  sigma {s}: GD/lin {mg/(vg*s*s):.3f} SMM/lin {ms/(vs*s*s):.3f} ratio {ms/mg:.4f}')
