# Checks for Chapter 6: exact two step identity for general f, third order defect, period two Chebyshev schedule.
import numpy as np
rng=np.random.default_rng(0); beta=np.log(1.5)/np.log(2)
d=5; X=rng.standard_normal((50,d)); yv=rng.integers(0,2,50)*2-1
def grad(w): return -(X.T@(yv/(1+np.exp(yv*(X@w)))))/50+0.1*w
L=np.linalg.eigvalsh(X.T@X/50).max()/4+0.1; al=1/L
x=rng.standard_normal(d); xp=x.copy(); xs=[x.copy()]
for t in range(40):
    if t%2==1: y=x+beta*(x-xp); new=y-al*grad(y)
    else: new=x-al*grad(x)
    xp,x=x,new; xs.append(x.copy())
z=xs[0].copy(); errs=[]; defect=[]
for j in range(20):
    u=z-(1+beta)*al*grad(z); z=u-al*grad(u); errs.append(np.linalg.norm(z-xs[2*j+2]))
    xt=xs[2*j]; x1=xs[2*j+1]; defect.append(np.linalg.norm(xs[2*j+2]-(x1-(1+beta)*al*grad(x1))))
print("max |even iterate - cyclic GD (1+b)a then a| :",max(errs))
print("defect of reversed order, first steps:",defect[:4])
# period two Chebyshev on [mu,L]
for kappa in [10,100,1000]:
    mu,LL=1.0,float(kappa); lp=(LL+mu)/2+(LL-mu)/(2*np.sqrt(2)); lm=(LL+mu)/2-(LL-mu)/(2*np.sqrt(2))
    lam=np.linspace(mu,LL,200001); P=np.abs((1-lam/lp)*(1-lam/lm)).max()
    s=(kappa+1)/(kappa-1); print(kappa,"ratio of steps",lp/lm,"max|P|",P,"1/T2",1/(2*s*s-1),"per step",np.sqrt(P),
      "SMM per step at alpha=1/L", np.sqrt(np.abs((1-lam/LL)*(1-(1+beta)*lam/LL)).max()))
