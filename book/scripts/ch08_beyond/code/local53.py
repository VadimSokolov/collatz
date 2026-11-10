import numpy as np, math
from scipy.optimize import brentq
beta=math.log(1.5)/math.log(2); mu,L,c=0.05,1.0,1.5; alpha=1/L
g=lambda x: mu*x+(L-mu)*np.tanh(x-c)
xs=brentq(lambda x: mu*x+(L-mu)*math.tanh(x-c),-10,10); H=mu+(L-mu)/math.cosh(xs-c)**2; a=alpha*H
Vb=beta+2/(2+beta)
def Vgd(a,s): return alpha**2*s**2/(a*(2-a))
def Vsmm(a,s):
    rho=(1-a)*(1-(1+beta)*a); ve=alpha**2*s**2*((1-a)**2*(1+beta)**2+1)/(1-rho**2); vo=(1-a)**2*ve+alpha**2*s**2
    return ve,vo,0.5*(ve+vo)
print('x*',xs,'H',H,'a',a,'V(beta)',Vb,'kappa',L/mu)
rng=np.random.default_rng(0); nch=4000; T=6000; burn=2000
for s in [1.0,0.5,0.2,0.05]:
    # GD
    x=np.full(nch,xs); acc=0;n=0
    for t in range(T):
        x=x-alpha*(g(x)+s*rng.standard_normal(nch))
        if t>=burn: acc+=((x-xs)**2).mean(); n+=1
    mgd=acc/n
    # SMM det
    x=np.full(nch,xs); xp=x.copy(); ae=ao=0; n=0
    for t in range(T):
        if t%2==1:
            y=x+beta*(x-xp); xn=y-alpha*(g(y)+s*rng.standard_normal(nch))
        else:
            xn=x-alpha*(g(x)+s*rng.standard_normal(nch))
        xp,x=x,xn
        if t>=burn:
            if (t+1)%2==0: ae+=((x-xs)**2).mean()
            else: ao+=((x-xs)**2).mean()
            n+=1
    me,mo=ae/(n/2),ao/(n/2); ms=0.5*(me+mo)
    ve,vo,vs=Vsmm(a,s); eps=math.sqrt(alpha)*s*(L-mu)/mu**1.5
    # worst-case bound (class): a=alpha*mu
    wb=Vsmm(alpha*mu,s)[2]
    print(f'sigma {s}: GD mse/lin {mgd/Vgd(a,s):.3f}  SMM mse/lin {ms/vs:.3f}  ratio SMM/GD {ms/mgd:.4f} (lin {vs/Vgd(a,s):.4f}, V(beta) {Vb:.4f})  SMM mse / class bound {ms/wb:.3f}')
# check worst-case ratio bound numerically over c1 for small kappa regime
import itertools
worst=0
for b in np.linspace(0,0.999,200):
    for am in np.linspace(0.001,1,400):   # am=alpha*mu, with alpha L =1
        c1=1-am; c2=max(1-(1+b)*am, b)
        q=c1*c2; ve=((1+b)**2*c1**2+1)/(1-q*q); vo=c1*c1*ve+1; r=0.5*(ve+vo)*(am*(2-am))
        worst=max(worst,r/(b+2/(2+b)))
print('max over (beta,alpha mu) of class bound /(V(beta) V_GD(alpha mu)) at alpha=1/L:',worst)
